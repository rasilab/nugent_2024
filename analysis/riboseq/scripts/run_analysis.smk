# useful libraries
import os
import pandas as pd
import re
import itertools as it


# configuration specific to this analysis
sample_annotations = (
  pd.read_table("../annotations/sample_annotations.csv", sep=",", dtype=object)
  .set_index('sample_id')
)

sample_ids = sample_annotations.index.tolist()


rule all:
    input:
        trimmed = expand('../data/trim/{sample_id}.fastq', sample_id=sample_ids),
        norrna = expand('../data/norrna/{sample_id}.fastq', sample_id=sample_ids),
        alignments = expand('../data/alignments/{sample_id}.bam', sample_id=sample_ids),
        plasmid_alignments = expand('../data/alignments/{sample_id}_plasmid.bam', sample_id=sample_ids),
        tx_counts = expand('../data/tx_counts/{sample_id}.tsv', sample_id=sample_ids),


def get_fastq_files_from_srr(srr):
  """This function returns the names of fastq files for parallel-fastq-dump
  """
  sample_id = sample_annotations.loc[sample_annotations['srr'] == srr].index[0]
  n_reads = 1
  filenames = [f'../data/fastq/{srr}_{read}.fastq' for read in range(1,n_reads+1)]
  return filenames


def get_split_read_files_input(wildcards):
  """This function returns the names of fastq files for combining them
  """
  srr = sample_annotations.loc[wildcards.sample_id, 'srr']
  # Get the results from the checkpoint, this is done here simply to make
  # the function depend on the checkpoint, but is not used in the function
  checkpoint_output = checkpoints.get_fastq.get(srr=srr).output[0]
  n_reads = 1
  filenames = [f'../data/fastq/{srr}_{read}.fastq' for read in range(1,n_reads+1)]
  return filenames


rule download_sra:
  """Download SRA file from NCBI database"""
  input:
  output:
    '../data/sra/{srr}.sra'
  threads: 1 
  container: "docker://ghcr.io/rasilab/sratools:3.0.8"
  shell:
    """
    set +e # continue if there is an error code
    prefetch {wildcards.srr} --output-file {output} 
    exitcode=$?
    if [ $exitcode -eq 1 ]
    then
      exit 1
    else
      exit 0
    fi    
    """


checkpoint get_fastq:
  """Convert SRA to variable fastq"""
  input: '../data/sra/{srr}.sra'
  output: '../data/fastq/{srr}_1.fastq'
  params:
    directory = '../data/fastq'
  threads: 36
  container: "docker://ghcr.io/rasilab/parallel_fastq_dump:0.6.7"
  shell:
    """
    set +e # continue if there is an error code
    parallel-fastq-dump \
      --sra-id {input} \
      --threads {threads} \
      --outdir {params.directory} \
      --tmpdir {params.directory} \
      --split-files
    exitcode=$?
    if [ $exitcode -eq 1 ]
    then
      exit 1
    else
      exit 0
    fi    
    """


rule trim_linker:
  """Remove extra sequences before aligning"""
  input: get_split_read_files_input
  params:
    adapter = lambda wildcards, output: sample_annotations.loc[wildcards.sample_id, 'adapter'],
    trim5 = lambda wildcards, output: sample_annotations.loc[wildcards.sample_id, 'trim5'],
    trim3 = lambda wildcards, output: sample_annotations.loc[wildcards.sample_id, 'trim3'],
    minimum_length = 22
  output: '../data/trim/{sample_id}.fastq'
  threads: 36
  log:
    '../data/trim/{sample_id}.log'
  container: "docker://ghcr.io/rasilab/cutadapt:4.4"
  shell:
      """
      cutadapt \
      --cores={threads} \
      --adapter={params.adapter} \
      --cut={params.trim5} \
      --cut=-{params.trim3} \
      --minimum-length={params.minimum_length} \
      --match-read-wildcards \
      --output {output} \
      {input} \
      &> {log}
      """


rule download_ribosomal_rrna_sequences:
  """Download ribosomal RNA sequences"""
  output:
    rrna_28s = '../data/rrna/28s.rrna.fa',
    rrna_18s = '../data/rrna/18s.rrna.fa',
    rrna_5s = '../data/rrna/5s.rrna.fa',
    rrna_58s = '../data/rrna/5.8s.rrna.fa'
  params:
    rrna_28s = "NR_003287.2",
    rrna_18s = "NR_003286.3",
    rrna_5s = "NR_023363.1",
    rrna_58s = "NR_003285.2"
  container: "docker://ghcr.io/rasilab/entrez-direct:16.2"
  shell:
    """
    esearch -db nucleotide -query {params.rrna_28s} | efetch -format fasta > {output.rrna_28s}
    esearch -db nucleotide -query {params.rrna_18s} | efetch -format fasta > {output.rrna_18s}
    esearch -db nucleotide -query {params.rrna_5s} | efetch -format fasta > {output.rrna_5s}
    esearch -db nucleotide -query {params.rrna_58s} | efetch -format fasta > {output.rrna_58s}
    """


rule make_bowtie_rrna_index:
  """Make Bowtie2 index for ribosomal RNA"""
  input:
    rrna_28s = '../data/rrna/28s.rrna.fa',
    rrna_18s = '../data/rrna/18s.rrna.fa',
    rrna_5s = '../data/rrna/5s.rrna.fa',
    rrna_58s = '../data/rrna/5.8s.rrna.fa'
  output:
    index = '../data/bowtie/rrna.1.ebwt'
  threads: 12
  params:
    index = '../data/bowtie/rrna'
  container: "docker://ghcr.io/rasilab/bowtie:1.3.1"
  shell:
    """
    bowtie-build \
    --threads {threads} \
    {input.rrna_28s},{input.rrna_18s},{input.rrna_5s},{input.rrna_58s} \
    {params.index} &> {params.index}.log
    """


rule remove_rrna:
  """Remove contaminant reads aligning to ribosomal rRNA"""
  input:
    reads = '../data/trim/{sample_id}.fastq',
    index = '../data/bowtie/rrna.1.ebwt'
  params:
    rrna_index = '../data/bowtie/rrna'
  output:
    '../data/norrna/{sample_id}.fastq'
  log:
    '../data/norrna/{sample_id}.log'
  threads: 36
  container: "docker://ghcr.io/rasilab/bowtie:1.3.1"
  shell:
    """
    bowtie \
    --threads {threads} \
    --un {output} \
    -x {params.rrna_index} \
    {input.reads} \
    1> /dev/null \
    2> {log} 
    """


rule download_mane_tx_fasta:
    input:
    output: "../data/mane/MANE.GRCh38.v1.3.ensembl_rna.fna"
    singularity: "docker://ghcr.io/rasilab/python:1.0.0"
    params:
      annotations_url = 'https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.3/MANE.GRCh38.v1.3.ensembl_rna.fna.gz'
    shell:
        """
        wget -O - {params.annotations_url} | gunzip -c 1> {output} 2> {output}.log
        """


rule make_transcriptome_bowtie_index:
  input: "../data/mane/MANE.GRCh38.v1.3.ensembl_rna.fna",
  output: "../data/mane/MANE.GRCh38.v1.3.ensembl_rna.1.ebwt",
  params: 
    index = "../data/mane/MANE.GRCh38.v1.3.ensembl_rna",
  threads: 36
  container: "docker://ghcr.io/rasilab/bowtie:1.3.1"
  shell:
    """
    bowtie-build \
    --threads {threads} \
    {input} \
    {params.index} &> {params.index}.log
    """


rule align_transcriptome:
  input:
    reads = "../data/norrna/{sample}.fastq",
    index = "../data/mane/MANE.GRCh38.v1.3.ensembl_rna.1.ebwt"
  output:
    aligned = "../data/alignments/{sample}.sam"
  log: "../data/alignments/{sample}.bowtie.log"
  params:
    index =  "../data/mane/MANE.GRCh38.v1.3.ensembl_rna"
  threads: 36
  container: "docker://ghcr.io/rasilab/bowtie:1.3.1"
  shell:
    """
    bowtie \
    --norc \
    --sam \
    --no-unal \
    --threads {threads} \
    -x {params.index} \
    {input.reads} \
    1> {output.aligned} 2> {log}
    """


rule make_plasmid_bowtie_index:
  input: "../annotations/plasmids/pHPHS232_pHPHS800_pPNHS189.fa",
  output: "../data/plasmid_bowtie/pHPHS232_pHPHS800_pPNHS189.1.ebwt",
  params: 
    index = "../data/plasmid_bowtie/pHPHS232_pHPHS800_pPNHS189",
  threads: 36
  container: "docker://ghcr.io/rasilab/bowtie:1.3.1"
  shell:
    """
    bowtie-build \
    --threads {threads} \
    {input} \
    {params.index} &> {params.index}.log
    """


rule align_plasmids:
  input:
    reads = "../data/norrna/{sample}.fastq",
    index = "../data/plasmid_bowtie/pHPHS232_pHPHS800_pPNHS189.1.ebwt"
  output:
    aligned = "../data/alignments/{sample}_plasmid.sam"
  log: "../data/alignments/{sample}_plasmid.bowtie.log"
  params:
    index =  "../data/plasmid_bowtie/pHPHS232_pHPHS800_pPNHS189"
  threads: 36
  container: "docker://ghcr.io/rasilab/bowtie:1.3.1"
  shell:
    """
    bowtie \
    --norc \
    --sam \
    --no-unal \
    --threads {threads} \
    -x {params.index} \
    {input.reads} \
    1> {output.aligned} 2> {log}
    """


rule sort_and_index_alignments:
    """Sort and index alignments
    """
    input:
        sam = '../data/alignments/{sample}.sam'
    output:
        bam = '../data/alignments/{sample}.bam',
        bai = '../data/alignments/{sample}.bam.bai',
    threads: 36
    singularity: "docker://ghcr.io/rasilab/samtools:1.16.1"
    shell:
        """
        samtools view -@ {threads} -b {input.sam} > {output.bam}.unsorted
        samtools sort -@ {threads} {output.bam}.unsorted > {output.bam}
        samtools index -@ {threads} {output.bam}
        rm {output.bam}.unsorted
        """


rule count_transcript_reads:
    """Count reads aligning to each transcript
    """
    input:
        bam = '../data/alignments/{sample}.bam',
    output:
        counts = '../data/tx_counts/{sample}.tsv',
    threads: 1
    container: "docker://ghcr.io/rasilab/samtools:1.16.1"
    shell:
        """
        samtools idxstats {input.bam} | cut -f 1,3 | sort -k2nr > {output.counts}
        """
