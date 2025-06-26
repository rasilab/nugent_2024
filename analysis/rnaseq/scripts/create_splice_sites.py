#!/usr/bin/env python3

# Adapted from HISAT2 splice site extraction script
# Creates annotated splice sites and comprehensive annotation file

from sys import stderr, exit, argv
from collections import defaultdict as dd

def extract_splice_sites(gtf_file, output_prefix):
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
            genes[gene_id].append(transcript_id)
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

    # Calculate annotated junctions (only consecutive exons)
    junctions = []
    for tran_id, [chrom, strand, exons, gene_id, gene_name] in trans.items():
        for i in range(len(exons)-1):
            # Only consecutive exons (annotated splice junctions)
            junction = {
                'seqnames': chrom,
                'start': exons[i][1] + 1,  # End of exon i + 1
                'end': exons[i+1][0] - 1,  # Start of exon i+1 - 1
                'strand': strand,
                'transcript_id': tran_id,
                'gene_id': gene_id,
                'gene_name': gene_name,
                'junction_type': 'annotated',
                'n_exon_skip': 0
            }
            junctions.append(junction)

    # Write comprehensive annotation file (TSV)
    import pandas as pd
    sj_df = pd.DataFrame(junctions)
    sj_df.to_csv(f'{output_prefix}.tsv', sep='\t', index=False)
    
    # Write simple splice sites file for STAR
    with open(f'{output_prefix}.ss', 'w') as f:
        for _, row in sj_df.iterrows():
            f.write(f"{row['seqnames']}\t{row['start']}\t{row['end']}\t{row['strand']}\n")
    
    print(f"Created {len(junctions)} annotated splice junctions from {len(trans)} transcripts", file=stderr)
    return len(junctions)

if __name__ == '__main__':
    if len(argv) != 3:
        print("Usage: python create_splice_sites.py input.gtf output_prefix", file=stderr)
        exit(1)
    
    gtf_file = argv[1]
    output_prefix = argv[2]
    
    with open(gtf_file, 'r') as f:
        extract_splice_sites(f, output_prefix)