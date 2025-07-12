"""Workflow for linking sgRNAs to barcodes

  :Author: Arvind Rasi Subramaniam (rasi@fredhutch.org)
"""

# useful libraries
import os
import pandas as pd
import re
import itertools as it

# Container directory path
container_dir = config.get('container_dir', '../../../.env/singularity_cache')


# configuration specific to this analysis
sample_annotations = pd.read_table("../annotations/sample_annotations.csv", 
                                   sep=",", comment = "#", dtype=object).set_index('sample_id')
print(sample_annotations)
read_count_cutoff =  5


# these rules are run locally
localrules: all

# Rules ----------------------------------------------------------------------

rule all:
  """List of all files we want at the end
  """
  input:
    sra = expand("../data/sra/{srr}.sra", srr=sample_annotations['srr']),
    fastq = expand("../data/fastq/{srr}_{read}.fastq", srr=sample_annotations['srr'], read=range(1,4)),
    insert_barcode_counts = expand('../data/insert_barcode_counts/{sample_id}.csv', 
      sample_id=sample_annotations.index),
    annotated_insert_barcode_counts = expand('../data/annotated_insert_barcode_counts/{sample_id}.csv', 
      sample_id=sample_annotations.index),
    ref_vs_ref_align = expand('../data/ref_vs_ref_alignments/{sample_id}/alignment_barcode{barcode_num}.bam',
      sample_id=sample_annotations.index, barcode_num=[1,2]),
    filtered_barcodes = expand('../data/filtered_barcodes/{sample_id}.csv',
      sample_id=sample_annotations.index),
   

def get_split_read_files_input(wildcards):
  """This function returns the names of R1,R2,R3 files for combining them
  """
  srr = sample_annotations.loc[wildcards.sample_id, 'srr']
  filenames = [f'../data/fastq/{srr}_{read}.fastq' for read in range(1,4)]
  return filenames


rule download_sra:
  """Download SRA file from NCBI database"""
  input:
  output:
    '../data/sra/{srr}.sra'
  threads: 1 
  singularity: f"{container_dir}/sratools_3.0.8.sif"
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


rule get_fastq:
  """Download fastq from SRA"""
  input: '../data/sra/{srr}.sra'
  output:
        ['../data/fastq/{srr}_' + f'{read}.fastq' for read in range(1,4)]
  params:
    directory = '../data/fastq'
  threads: 36
  singularity: f"{container_dir}/parallel_fastq_dump_0.6.7.sif"
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


rule count_insert_barcodes:
  """Extract and count distinct insert-barcode combinations
  """
  input: get_split_read_files_input
  output: '../data/insert_barcode_counts/{sample_id}.csv'
  params:
    insert1_read = lambda w: sample_annotations.loc[w.sample_id, 'insert1_read'],
    insert1_start = lambda w: sample_annotations.loc[w.sample_id, 'insert1_start'],
    insert1_length = lambda w: sample_annotations.loc[w.sample_id, 'insert1_length'],
    insert2_read = lambda w: sample_annotations.loc[w.sample_id, 'insert2_read'],
    insert2_start = lambda w: sample_annotations.loc[w.sample_id, 'insert2_start'],
    insert2_length = lambda w: sample_annotations.loc[w.sample_id, 'insert2_length'],
    barcode1_read = lambda w: sample_annotations.loc[w.sample_id, 'barcode1_read'],
    barcode1_start = lambda w: sample_annotations.loc[w.sample_id, 'barcode1_start'],
    barcode1_length = lambda w: sample_annotations.loc[w.sample_id, 'barcode1_length'],
    barcode2_read = lambda w: sample_annotations.loc[w.sample_id, 'barcode2_read'],
    barcode2_start = lambda w: sample_annotations.loc[w.sample_id, 'barcode2_start'],
    barcode2_length = lambda w: sample_annotations.loc[w.sample_id, 'barcode2_length'],
  log: '../data/insert_barcode_counts/{sample_id}.log'
  singularity: f"{container_dir}/python_1.0.0.sif"
  shell: 
    """
    set +o pipefail;
    export TMPDIR=$(pwd);
    paste {input} | 
    awk '
      BEGIN {{
        OFS = ",";
        SUBSEP = ",";
        print "counts","insert1","insert2","barcode1","barcode2"
      }}
      NR%4 == 2  {{ 
        counts[\
              substr(${params.insert1_read}, {params.insert1_start}, {params.insert1_length}), \
              substr(${params.insert2_read}, {params.insert2_start}, {params.insert2_length}), \
              substr(${params.barcode1_read}, {params.barcode1_start}, {params.barcode1_length}), \
              substr(${params.barcode2_read}, {params.barcode2_start}, {params.barcode2_length}) \
              ]++
      }}
      END {{
        for (key in counts) print counts[key],key | "sort -k1nr"
      }}
    ' 1> {output}  2> {log}
    """


rule subset_to_annotated_inserts:
  """Subset insert-barcode counts to only annotated inserts
  """
  input:
    count_file = '../data/insert_barcode_counts/{sample_id}.csv',
    insert_annotations = '../annotations/insert_annotations.csv',
  output:
    annotated_insert_count_file = '../data/annotated_insert_barcode_counts/{sample_id}.csv'
  params:
    # column in insert annotations file that contains the insert number and sequence
    insert_num_column = 1,
    insert1_seq_column = 2,
    insert2_seq_column = 3,
  log: '../data/annotated_insert_barcode_counts/{sample_id}.log'
  singularity: f"{container_dir}/python_1.0.0.sif"
  shell:
    """
    awk '
      # set up parameters and header
      BEGIN {{
        FS=",";
        OFS=",";
        SUBSEP=",";
        barcode_num=1;
        insert_num_col={params.insert_num_column}; 
        insert1_seq_col={params.insert1_seq_column};
        insert2_seq_col={params.insert2_seq_column};
        print "count","barcode_num","insert_num","barcode1","barcode2"
      }}
      # read in all insert annotations
      NR == FNR {{
        inserts[$insert1_seq_col,$insert2_seq_col] = $insert_num_col; 
        next
      }}
      # write if insert is present
      {{
        if ($2 SUBSEP $3 in inserts) {{
          print $1,barcode_num,inserts[$2,$3],$4,$5;
          barcode_num++
        }}
      }}
      ' \
      {input.insert_annotations} {input.count_file} 1> {output.annotated_insert_count_file} 2> {log}
    """


rule align_barcodes1_against_themselves:
  """Align barcodes against themselves to find multialigners
  """
  input:
    '../data/annotated_insert_barcode_counts/{sample_id}.csv'
  output:
    sam = temp('../data/ref_vs_ref_alignments/{sample_id}/alignment_barcode1.sam'),
    bam = '../data/ref_vs_ref_alignments/{sample_id}/alignment_barcode1.bam',
    fasta = '../data/ref_vs_ref_alignments/{sample_id}/reference_barcode1.fasta',
  log:
    align = '../data/ref_vs_ref_alignments/{sample_id}/align_barcode1.log',
    build = '../data/ref_vs_ref_alignments/{sample_id}/build_barcode1.log',
  params:
    bowtie_index = '../data/ref_vs_ref_alignments/{sample_id}/reference_barcode1',
    barcode_count_column = 1,
    barcode_num_column = 2,
    barcode_seq_column = 4,
    read_count_cutoff = read_count_cutoff,
  threads: 36
  singularity: f"{container_dir}/bowtie2_2.4.5.sif"
  shell:
    """
    # write the input file to a fasta file of barcodes with name as barcode_num col from input_file
    awk -F, 'NR > 1 {{if (${params.barcode_count_column} >= {params.read_count_cutoff}) print ">" ${params.barcode_num_column} "\\n" ${params.barcode_seq_column}}}' {input} > {output.fasta}
    # create a bowtie reference of the barcodes
    bowtie2-build {output.fasta} {params.bowtie_index} 2> {log.build}
    # align against itself
    bowtie2 --threads {threads} -L 19 -N 1 --all --norc --no-unal -f -x {params.bowtie_index} -U {output.fasta} > {output.sam}  2> {log.align}
    # convert to BAM
    samtools view -@ {threads} -b {output.sam} > {output.bam}.tmp
    # sort
    samtools sort -@ {threads} -m 20M {output.bam}.tmp > {output.bam}
    sleep 10
    # index
    samtools index -@ {threads} {output.bam}
    # remove unsorted bam
    rm {output.bam}.tmp
    """


rule align_barcodes2_against_themselves:
  """Align barcodes against themselves to find multialigners
  """
  input:
    '../data/annotated_insert_barcode_counts/{sample_id}.csv'
  output:
    sam = temp('../data/ref_vs_ref_alignments/{sample_id}/alignment_barcode2.sam'),
    bam = '../data/ref_vs_ref_alignments/{sample_id}/alignment_barcode2.bam',
    fasta = '../data/ref_vs_ref_alignments/{sample_id}/reference_barcode2.fasta',
  log:
    align = '../data/ref_vs_ref_alignments/{sample_id}/align_barcode2.log',
    build = '../data/ref_vs_ref_alignments/{sample_id}/build_barcode2.log',
  params:
    bowtie_index = '../data/ref_vs_ref_alignments/{sample_id}/reference_barcode2',
    barcode_count_column = 1,
    barcode_num_column = 2,
    barcode_seq_column = 5,
    read_count_cutoff = read_count_cutoff,
  threads: 36
  singularity: f"{container_dir}/bowtie2_2.4.5.sif"
  shell:
    """
    # write the input file to a fasta file of barcodes with name as barcode_num col from input_file
    awk -F, 'NR > 1 {{if (${params.barcode_count_column} >= {params.read_count_cutoff}) print ">" ${params.barcode_num_column} "\\n" ${params.barcode_seq_column}}}' {input} > {output.fasta}
    # create a bowtie reference of the barcodes
    bowtie2-build {output.fasta} {params.bowtie_index} 2> {log.build}
    # align against itself
    bowtie2 --threads {threads} -L 19 -N 1 --all --norc --no-unal -f -x {params.bowtie_index} -U {output.fasta} > {output.sam}  2> {log.align}
    # convert to BAM
    samtools view -@ {threads} -b {output.sam} 1> {output.bam}.tmp 2> {log.align}
    # sort
    samtools sort -@ {threads} -m 20M {output.bam}.tmp 1> {output.bam} 2> {log.align}
    sleep 10
    # index
    samtools index -@ {threads} {output.bam}
    # remove unsorted bam
    rm {output.bam}.tmp
    """


rule filter_barcodes:
  """Filter barcodes to remove clashes and sequencing errors and produce a final list
  """
  input:
    bam1 = '../data/ref_vs_ref_alignments/{sample_id}/alignment_barcode1.bam',
    bam2 = '../data/ref_vs_ref_alignments/{sample_id}/alignment_barcode2.bam',
    counts = '../data/annotated_insert_barcode_counts/{sample_id}.csv',
    notebook = 'filter_barcodes.ipynb',
  output:
    '../data/filtered_barcodes/{sample_id}.csv',
  params:
    read_count_cutoff = read_count_cutoff
  log:
    '../data/filtered_barcodes/{sample_id}.log',
  singularity: f"{container_dir}/r_1.0.0.sif"
  shell:
    """
    jupyter nbconvert --to script --ExecutePreprocessor.kernel_name=ir {input.notebook}

    notebook={input.notebook}
    script="${{notebook/.ipynb/.r}}"

    Rscript ${{script}} {input.bam1} {input.bam2} {input.counts} {params.read_count_cutoff} {output} &> {log}
    """