"""
Snakemake workflow for RNA-seq alignment using STAR aligner

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


def get_fastq_files(wildcards):
    """Get input files based on sample_id
    """
    fastq = [f'../data/fastq/{filename}' for filename in filter(lambda x: re.search(f'{wildcards.sample_id}_', x) and re.search('_R[12]_', x), sorted(os.listdir('../data/fastq/')))]
    return fastq


rule all:
    input:
        genome = [f"../data/ensembl/genome/Homo_sapiens.GRCh38.dna_rm.chromosome.{chr}.fa"
                    for chr in list(range(1,23)) + ['X', 'MT']],
        annotations = "../data/ensembl/Homo_sapiens.GRCh38.108.gtf",
        genome_splice_sites = "../data/ensembl/Homo_sapiens.GRCh38.108.all.ss",
        genome_splice_annotations = "../data/ensembl/Homo_sapiens.GRCh38.108.all.ss.tsv",
        genome_index = "../data/star/Homo_sapiens.GRCh38genome.108tran/genomeParameters.txt",
        genome_alignments = [f'../data/alignments/genome_{sample_id}/Aligned.out.bam' for sample_id in sample_ids],
        genome_intron_counts = [f'../data/intron_counts/genome_{sample_id}.csv' for sample_id in sample_ids],


rule download_genome:
    """Download genome from ensembl"""
    input:
    output: "../data/ensembl/genome/Homo_sapiens.GRCh38.dna_rm.chromosome.{chr}.fa"
    log: "../data/ensembl/genome/wget.{chr}.log"
    singularity: "docker://ghcr.io/rasilab/python:1.0.0"
    params:
        filename = "Homo_sapiens.GRCh38.dna_rm.chromosome.{chr}.fa.gz",
        genome_url = "https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/"
    shell:
        """
        wget -O - {params.genome_url}{params.filename} | gunzip -c 1> {output} 2> {log}
        """


rule download_ensembl_annotations:
    """Download ensembl annotations"""
    input:
    output: "../data/ensembl/Homo_sapiens.GRCh38.108.gtf",
    log:
        ensembl_log = "../data/ensembl/wget.log"
    singularity: "docker://ghcr.io/rasilab/python:1.0.0"
    params:
        chr = list(range(1,23)) + ['X', 'MT'],
        annotations_url = "https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz"
    shell:
        """
        # use awk to extract only annotations in primary chromosomes
        wget -O - {params.annotations_url} | gunzip -c | 
        awk 'BEGIN {{split("{params.chr}", chr, " ")}} \
                   {{for (a in chr) if ($1 == chr[a]) print $0}}' \
                   1> {output} 2> {log}
        """


rule create_ensembl_splice_site_database:
    """Create ensembl splice site database for STAR
    """
    input:        
        annotations = "../data/ensembl/Homo_sapiens.GRCh38.108.gtf"
    output:
        annotated_splice_sites = "../data/ensembl/Homo_sapiens.GRCh38.108.annotated.ss",
        all_splice_sites = "../data/ensembl/Homo_sapiens.GRCh38.108.all.ss",
        splice_annotations = "../data/ensembl/Homo_sapiens.GRCh38.108.all.ss.tsv"
    singularity: "docker://ghcr.io/rasilab/python:1.0.0"
    shell:
        """
        python extract_splice_site_annotations.py {input.annotations} 1> {output.splice_annotations}
        awk -v OFS="\\t" 'NR > 1 {{print $1,$2,$3,$4}}' {output.splice_annotations} 1> {output.all_splice_sites}
        awk -v OFS="\\t" '$8 == "annotated" {{print $1,$2,$3,$4}}' {output.splice_annotations} 1> {output.annotated_splice_sites}
        """


rule create_genome_index_for_alignment:
    """Create genome index for alignment
    """
    input:
        genome = [f"../data/ensembl/genome/Homo_sapiens.GRCh38.dna_rm.chromosome.{chr}.fa"
                    for chr in list(range(1,23)) + ['X', 'MT']],
        gtf = "../data/ensembl/Homo_sapiens.GRCh38.108.gtf",
    output:
        index_file = "../data/star/Homo_sapiens.GRCh38genome.108tran/genomeParameters.txt"
    params:
        genomeDir = "../data/star/Homo_sapiens.GRCh38genome.108tran/"
    # singularity: "docker://ghcr.io/rasilab/star:2.7.11a"
    envmodules: 'STAR/2.7.10b-GCC-12.2.0'
    threads: 36
    shell:
        """
        STAR \
        --runThreadN {threads} \
        --runMode genomeGenerate \
        --sjdbGTFfile {input.gtf} \
        --limitSjdbInsertNsj 3000000 \
        --genomeFastaFiles {input.genome} \
        --genomeDir {params.genomeDir}
        """


rule genome_annotated_align:
    """Align reads against genome with annotated splice sites
    """
    input: 
        fastq = get_fastq_files,
        index_file = '../data/star/Homo_sapiens.GRCh38genome.108tran/genomeParameters.txt'
    output: 
        alignments = "../data/alignments/genome_{sample_id}/Aligned.out.sam",
        fastq_read1 = '../data/alignments/genome_{sample_id}/Unmapped.out.mate1', 
        fastq_read2 = '../data/alignments/genome_{sample_id}/Unmapped.out.mate2', 
    threads: 36
    params:
        genomeDir = "../data/star/Homo_sapiens.GRCh38genome.108tran",
        input_file_string = get_fastq_files,
        outFileNamePrefix = "../data/alignments/genome_{sample_id}/"
    envmodules: 'STAR/2.7.10b-GCC-12.2.0'
    shell:
        """
        STAR \
        --runThreadN {threads} \
        --runMode alignReads \
        --alignSJoverhangMin 300 \
        --alignSJDBoverhangMin 6 \
        --outSAMmultNmax 1 \
        --quantMode GeneCounts \
        --outReadsUnmapped Fastx \
        --genomeDir  {params.genomeDir} \
        --readFilesIn {params.input_file_string} \
        --readFilesCommand zcat \
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


rule calculate_intron_counts_genome:
    input: 
        sj_annotations = "../data/ensembl/Homo_sapiens.GRCh38.108.all.ss.tsv",
        aln = "../data/alignments/genome_{sample_name}/Aligned.out.bam",
        notebook = 'get_intron_coverage.ipynb'
    output: '../data/intron_counts/genome_{sample_name}.csv'
    log: '../data/intron_counts/genome_{sample_name}.log'
    threads: 1
    singularity: "docker://ghcr.io/rasilab/r:1.0.0"
    shell:
        """
        export JUPYTER_DATA_DIR=$(pwd)
        export JUPYTER_CONFIG_DIR=$(pwd)
        jupyter nbconvert --to script --ExecutePreprocessor.kernel_name=ir {input.notebook}
        notebook={input.notebook}
        script="${{notebook/.ipynb/.r}}"
        Rscript ${{script}} \
            {input.sj_annotations} \
            {input.aln} \
            {output} \
            &> {log}
        """