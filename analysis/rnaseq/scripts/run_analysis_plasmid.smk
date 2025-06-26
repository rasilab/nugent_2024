"""
Snakemake workflow for analysis of RNA-seq data mapping to integrated plasmids

Arvind Rasi Subramaniam (rasi@fredhutch.org)
"""

import os
import pandas as pd
import itertools as it

# configuration specific to this analysis
sample_annotations = pd.read_table("../annotations/sample_annotations.csv", 
                                   sep=",", comment="#", dtype=object).set_index('sample_id')
print(sample_annotations)

sample_ids = [x for x in sample_annotations.index.tolist()]
study_list = sample_annotations['study'].unique()

plasmids = [x.split('.')[0] for x in os.listdir("../annotations/plasmids/") if x.endswith('.cleaned.gtf')]


rule all:
    input:
        plasmid_splice_sites = expand("../annotations/plasmids/{plasmid}.ss", plasmid=plasmids),
        plasmid_splice_annotations = expand("../annotations/plasmids/{plasmid}.ss.tsv", plasmid=plasmids),
        plasmid_index = expand("../data/star/{plasmid}/genomeParameters.txt", plasmid=plasmids),
        plasmid_alignments = [f'../data/alignments/plasmids_{sample_id}/Aligned.out.bam' for sample_id in sample_ids],
        plasmid_intron_counts = [f'../data/intron_counts/plasmids_{sample_id}.csv' for sample_id in sample_ids],
        plasmid_coverage_figures = ["../figures/globin_cvg.png"],


def get_split_read_files_input(wildcards):
    """This function returns the names of fastq files for combining them
    """
    srr = sample_annotations.loc[wildcards.sample_id, 'srr']
    n_reads = 2
    filenames = [f'../data/fastq/{srr}_{read}.fastq' for read in range(1, n_reads+1)]
    return filenames


rule create_plasmid_splice_site_database:
    """Create splice site database for STAR 
    """
    input:        
        annotations = "../annotations/plasmids/{plasmid}.cleaned.gtf"
    output:
        splice_sites = "../annotations/plasmids/{plasmid}.ss",
        splice_annotations = "../annotations/plasmids/{plasmid}.ss.tsv",
    singularity: "docker://ghcr.io/rasilab/python:1.0.0"
    shell:
        """
        python extract_all_possible_splice_sites.py {input.annotations} 1> {output.splice_sites}
        python extract_splice_site_annotations.py {input.annotations} 1> {output.splice_annotations}
        """


rule create_plasmid_index_for_alignment:
    """Create plasmid index for alignment
    """
    input:
        genome = "../annotations/plasmids/{plasmid}.fa",
        gtf = "../annotations/plasmids/{plasmid}.cleaned.gtf",
        splice_sites = "../annotations/plasmids/{plasmid}.ss",
    output:
        index_file = '../data/star/{plasmid}/genomeParameters.txt'
    params:
        genomeDir = "../data/star/{plasmid}/",
        genomeSAindexNbases = 6
    singularity: "docker://ghcr.io/rasilab/star:2.7.11a"
    threads: 36
    shell:
        """
        STAR \
            --runThreadN {threads} \
            --runMode genomeGenerate \
            --sjdbGTFfile {input.gtf} \
            --sjdbFileChrStartEnd {input.splice_sites} \
            --genomeSAindexNbases {params.genomeSAindexNbases} \
            --quantMode GeneCounts \
            --genomeFastaFiles {input.genome} \
            --genomeDir {params.genomeDir}
        """        


rule plasmid_align:
    """Align reads against plasmid
    """
    input: 
        fastq_read1 = get_split_read_files_input,
        index_file = lambda w: f"../data/star/{sample_annotations.loc[w.sample_id,'plasmid']}/genomeParameters.txt"
    output: 
        alignments = "../data/alignments/plasmids_{sample_id}/Aligned.out.sam",
    threads: 36
    params:
        genomeDir = lambda w: f"../data/star/{sample_annotations.loc[w.sample_id, 'plasmid']}/",
        input_file_string = get_split_read_files_input,
        outFileNamePrefix = "../data/alignments/plasmids_{sample_id}/"
    singularity: "docker://ghcr.io/rasilab/star:2.7.11a"
    shell:
        """
        STAR \
            --runThreadN {threads} \
            --runMode alignReads \
            --alignSJoverhangMin 300 \
            --alignSJDBoverhangMin 6 \
            --outSAMmultNmax 1 \
            --quantMode GeneCounts \
            --genomeDir  {params.genomeDir} \
            --readFilesIn {params.input_file_string} \
            --outFileNamePrefix {params.outFileNamePrefix}
        """


rule sort_and_index_alignments:
    """Sort and index alignments
    """
    input:
        sam = '../data/{folder}/{sample_name}.sam'
    output:
        bam = '../data/{folder}/{sample_name}.bam',
        bai = '../data/{folder}/{sample_name}.bam.bai',
    threads: 36
    singularity: "docker://ghcr.io/rasilab/samtools:1.16.1"
    shell:
        """
        # convert to BAM
        samtools view -@ {threads} -b {input.sam} > {output.bam}.tmp;
        # sort
        samtools sort -@ {threads} {output.bam}.tmp > {output.bam};
        # index
        samtools index -@ {threads} {output.bam}
        # remove unsorted bam
        rm {output.bam}.tmp
        """


rule calculate_intron_counts_plasmid:
    input: 
        sj_annotations = lambda w: f"../annotations/plasmids/{sample_annotations.loc[w.sample_name, 'plasmid']}.ss.tsv",
        aln = "../data/alignments/plasmids_{sample_name}/Aligned.out.bam",
        script = 'get_intron_coverage.R'
    output: '../data/intron_counts/plasmids_{sample_name}.csv'
    log: '../data/intron_counts/plasmids_{sample_name}.log'
    threads: 1
    singularity: "docker://ghcr.io/rasilab/r:1.0.0"
    shell:
        """
        Rscript {input.script} \
            {input.sj_annotations} \
            {input.aln} \
            {output} \
            &> {log}
        """


rule make_plasmid_rna_seq_coverage_plots:
    input:
        annotations = "../annotations/sample_annotations.csv",
        gtf = "../annotations/plasmids/pHPHS232_pHPHS800_pAS321.cleaned.gtf",
        alignments = [f'../data/alignments/plasmids_{sample_id}/Aligned.out.bam' for sample_id in sample_ids],
        script = "make_plasmid_rna_seq_coverage_plots.R"
    output:
        fig = "../figures/globin_cvg.png"
    log: "../logs/make_plasmid_rna_seq_coverage_plots.log"
    threads: 1
    singularity: "docker://ghcr.io/rasilab/r:1.0.0"
    shell:
        """
        Rscript {input.script} &> {log}
        """