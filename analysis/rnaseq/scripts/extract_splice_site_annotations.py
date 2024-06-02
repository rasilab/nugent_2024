#!/usr/bin/env python3

#
# Copyright 2015, Daehwan Kim <infphilo@gmail.com>
#
# This file is part of HISAT 2.
#
# HISAT 2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HISAT 2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HISAT 2.  If not, see <http://www.gnu.org/licenses/>.
#

from sys import stderr, exit
from collections import defaultdict as dd, Counter
from argparse import ArgumentParser, FileType


def extract_splice_sites(gtf_file, verbose=False):
    genes = dd(list)
    trans = {}

    # Parse valid exon lines from the GTF file into a dict by transcript_id
    for line in gtf_file:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if '#' in line:
            line = line.split('#')[0].strip()

        try:
            chrom, source, feature, left, right, score, \
                strand, frame, values = line.split('\t')
        except ValueError:
            continue
        left, right = int(left), int(right)

        if feature != 'exon' or left >= right:
            continue

        values_dict = {}
        for attr in values.split(';'):
            if attr:
                attr, _, val = attr.strip().partition(' ')
                values_dict[attr] = val.strip('"')

        if 'gene_id' not in values_dict or \
                'transcript_id' not in values_dict:
            continue

        transcript_id = values_dict['transcript_id']
        gene_id = values_dict['gene_id']
        gene_name = values_dict.get('gene_name', gene_id)
        if transcript_id not in trans:
            trans[transcript_id] = [chrom, strand, [[left, right]], gene_id, gene_name]
            genes[values_dict['gene_id']].append(transcript_id)
        else:
            trans[transcript_id][2].append([left, right])

    # Sort exons and merge where separating introns are <=5 bps
    for tran, [chrom, strand, exons, gene_id, gene_name] in trans.items():
            exons.sort()
            tmp_exons = [exons[0]]
            for i in range(1, len(exons)):
                if exons[i][0] - tmp_exons[-1][1] <= 5:
                    tmp_exons[-1][1] = exons[i][1]
                else:
                    tmp_exons.append(exons[i])
            trans[tran] = [chrom, strand, tmp_exons, gene_id, gene_name]

    # Calculate and print the unique junctions
    junctions = dict()
    for tran, [chrom, strand, exons, gene_id, gene_name] in trans.items():
        for i in range(0, len(exons)-1):
            for j in range(i + 1, len(exons)):
                if j == i + 1:
                    junction_type= 'annotated'
                else:
                    junction_type = 'skip'
                key = (chrom, exons[i][1] + 1, exons[j][0] - 1, strand)
                if ((key not in junctions) or (key in junctions and junctions[key][3] == 'skip' and junction_type == 'annotated')):
                    junctions[key] = (tran, gene_id, gene_name, junction_type, j - i - 1)
    print('\t'.join(
        ["seqnames", "start", "end", "strand", 
         "transcript_id", "gene_id", "gene_name", "junction_type", "n_exon_skip"]))
    for key in sorted(junctions):
        value = junctions[key]
        print('\t'.join([str(x) for x in key + value]))
        
    # Print some stats if asked
    if verbose:
        exon_lengths, intron_lengths, trans_lengths = \
            Counter(), Counter(), Counter()
        for chrom, strand, exons in trans.values():
            tran_len = 0
            for i, exon in enumerate(exons):
                exon_len = exon[1]-exon[0]+1
                exon_lengths[exon_len] += 1
                tran_len += exon_len
                if i == 0:
                    continue
                intron_lengths[exon[0] - exons[i-1][1]] += 1
            trans_lengths[tran_len] += 1

        print('genes: {}, genes with multiple isoforms: {}'.format(
                len(genes), sum(len(v) > 1 for v in genes.values())),
              file=stderr)
        print('transcripts: {}, transcript avg. length: {:.0f}'.format(
                len(trans), sum(trans_lengths.elements())//len(trans)),
              file=stderr)
        print('exons: {}, exon avg. length: {:.0f}'.format(
                sum(exon_lengths.values()),
                sum(exon_lengths.elements())//sum(exon_lengths.values())),
              file=stderr)
        print('introns: {}, intron avg. length: {:.0f}'.format(
                sum(intron_lengths.values()),
                sum(intron_lengths.elements())//sum(intron_lengths.values())),
              file=stderr)
        print('average number of exons per transcript: {:.0f}'.format(
                sum(exon_lengths.values())//len(trans)),
              file=stderr)


if __name__ == '__main__':
    parser = ArgumentParser(
        description='Extract splice junctions from a GTF file')
    parser.add_argument('gtf_file',
        nargs='?',
        type=FileType('r'),
        help='input GTF file (use "-" for stdin)')
    parser.add_argument('-v', '--verbose',
        dest='verbose',
        action='store_true',
        help='also print some statistics to stderr')

    args = parser.parse_args()
    if not args.gtf_file:
        parser.print_help()
        exit(1)
    extract_splice_sites(args.gtf_file, args.verbose)
