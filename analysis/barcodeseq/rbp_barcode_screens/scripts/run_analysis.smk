"""Workflow for counting barcodes

  :Author: Arvind Rasi Subramaniam
  :Date: 1 Jan 2024
"""

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

illumina_sample_annotations = (
  pd.read_table("../annotations/sample_annotations.csv", sep=",", dtype=object)
  .drop_duplicates(subset=['illumina_sample_id'])
  .set_index('illumina_sample_id')
)

# calculate reads we need for pre-processing
illumina_sample_annotations['n_reads'] = illumina_sample_annotations[['barcode1_read','subpool_barcode_read','umi_read']].max(axis=1).astype(int) 

# barcodes with UMI counts below this cutoff are discarded
umi_cutoff = 20

# mageck comparisons
mageck_comparisons = pd.read_table("../annotations/mageck_comparisons.csv",
    sep=",", comment="#", dtype=object)


# these rules are run locally
localrules: all

# Rules ----------------------------------------------------------------------

rule all:
  """List of all files we want at the end
  """
  input:
    raw_counts = expand('../data/barcode_umi_and_read_counts/{illumina_sample_id}.csv', 
      illumina_sample_id=illumina_sample_annotations.index),
    linked_counts = expand('../data/linked_barcode_counts/{illumina_sample_id}.csv', 
      illumina_sample_id=illumina_sample_annotations.index),
    subpool_counts = expand('../data/subpool_barcode_counts/{sample_id}.csv', 
      sample_id=sample_annotations.index),
    library_stats = expand('../data/library_statistics/{sample_id}.csv',
      sample_id=sample_annotations.index),
    insert_counts = expand('../data/insert_counts/{sample_id}.csv',
      sample_id=sample_annotations.index),
    insert_counts_random_partition = expand('../data/insert_counts_random_partition/{sample_id}.csv',
     sample_id=sample_annotations.index),
    mageck_inputs_random_partition = expand('../data/mageck_random_partition/{sample_name[0]}_vs_{sample_name[1]}/input_counts.tsv', 
                           sample_name=list(zip(mageck_comparisons['treatment'], mageck_comparisons['control']))),
    mageck_gene_summary_random_partition = expand('../data/mageck_random_partition/{sample_name[0]}_vs_{sample_name[1]}/mageck.gene_summary.tsv', 
                           sample_name=list(zip(mageck_comparisons['treatment'], mageck_comparisons['control']))),
    mageck_inputs = expand('../data/mageck/{sample_name[0]}_vs_{sample_name[1]}/input_counts.tsv', 
                           sample_name=list(zip(mageck_comparisons['treatment'], mageck_comparisons['control']))),
    mageck_gene_summary = expand('../data/mageck/{sample_name[0]}_vs_{sample_name[1]}/mageck.gene_summary.tsv', 
                           sample_name=list(zip(mageck_comparisons['treatment'], mageck_comparisons['control']))),
    mageck_sgrna_summary = expand('../data/mageck/{sample_name[0]}_vs_{sample_name[1]}/mageck.sgrna_summary.tsv', 
                           sample_name=list(zip(mageck_comparisons['treatment'], mageck_comparisons['control']))),
    gene_summary = '../data/mageck/gene_summary_table.csv.gz',
    sgrna_summary = '../data/mageck/sgrna_summary_table.csv.gz',


def get_fastq_files_from_srr(srr):
  """This function returns the names of fastq files for parallel-fastq-dump
  """
  illumina_sample_id = illumina_sample_annotations.loc[illumina_sample_annotations['srr'] == srr].index[0]
  n_reads = illumina_sample_annotations.loc[illumina_sample_id, 'n_reads']
  filenames = [f'../data/fastq/{srr}_{read}.fastq' for read in range(1,n_reads+1)]
  return filenames


def get_split_read_files_input(wildcards):
  """This function returns the names of fastq files for combining them
  """
  srr = illumina_sample_annotations.loc[wildcards.illumina_sample_id, 'srr']
  # Get the results from the checkpoint, this is done here simply to make
  # the function depend on the checkpoint, but is not used in the function
  checkpoint_output = checkpoints.get_fastq.get(srr=srr).output[0]
  n_reads = illumina_sample_annotations.loc[wildcards.illumina_sample_id, 'n_reads']
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


rule extract_barcodes_and_count_umi_read:
  """Extract barcodes and count UMIs and reads per barcode combination
  """
  input: get_split_read_files_input
  output: raw_count_file = '../data/barcode_umi_and_read_counts/{illumina_sample_id}.csv'
  params:
    barcode1_read = lambda w: illumina_sample_annotations.loc[w.illumina_sample_id, 'barcode1_read'],
    barcode1_start = lambda w: illumina_sample_annotations.loc[w.illumina_sample_id, 'barcode1_start'],
    barcode1_length = lambda w: illumina_sample_annotations.loc[w.illumina_sample_id, 'barcode1_length'],
    subpool_barcode_read = lambda w: illumina_sample_annotations.loc[w.illumina_sample_id, 'subpool_barcode_read'],
    subpool_barcode_start = lambda w: illumina_sample_annotations.loc[w.illumina_sample_id, 'subpool_barcode_start'],
    subpool_barcode_length = lambda w: illumina_sample_annotations.loc[w.illumina_sample_id, 'subpool_barcode_length'],
    umi_read = lambda w: illumina_sample_annotations.loc[w.illumina_sample_id, 'umi_read'],
    umi_start = lambda w: illumina_sample_annotations.loc[w.illumina_sample_id, 'umi_start'],
    umi_length = lambda w: illumina_sample_annotations.loc[w.illumina_sample_id, 'umi_length'],
  log: "../data/barcode_umi_and_read_counts/{illumina_sample_id}.log"
  container: 'docker://ghcr.io/rasilab/python:1.0.0'
  shell: 
    """
    set +o pipefail;
    export TMPDIR=$(pwd);
    # paste concatenates all input files line by line with tab separator
    paste {input} | awk '
    BEGIN {{ 
      OFS=",";
      SUBSEP=",";
      print "subpool_barcode", "barcode", "umi_count", "read_count"
    }}
    NR % 4 == 2 {{ 
      subpool_barcode = substr(${params.subpool_barcode_read}, {params.subpool_barcode_start}, {params.subpool_barcode_length});
      barcode = substr(${params.barcode1_read}, {params.barcode1_start}, {params.barcode1_length});
      umi = substr(${params.umi_read}, {params.umi_start}, {params.umi_length});
      umi_read_counts[subpool_barcode, barcode, umi]++;
    }}
    END {{
      for (e in umi_read_counts) {{
        split(e, a, SUBSEP);
        umi_counts[a[1],a[2]]++;
        read_counts[a[1],a[2]] += umi_read_counts[e]
      }}
      for (e in umi_counts) print e,umi_counts[e],read_counts[e] | "sort -t, -k3nr"
    }}
    ' 1> {output}  2> {log}
    """


rule subset_to_linked_barcodes:
  """Subset barcode-UMI counts to only those barcodes identified in linkage sequencing
  """
  input:
    raw_count_file = '../data/barcode_umi_and_read_counts/{illumina_sample_id}.csv',
    barcode_linkage_file = '../../rbp_sgrna_barcode_linkage/data/filtered_barcodes/182p1.csv'
  output:
    linked_count_file = '../data/linked_barcode_counts/{illumina_sample_id}.csv'
  params:
    barcode_seq_column = 2,
    linkage_ref_barcode_number = lambda w: illumina_sample_annotations.loc[w.illumina_sample_id, 'linkage_ref_barcode_number'],
    offset_to_linked_barcode = 2,
    linked_insert_num_column = 1,
    linked_barcode_num_column = 2,
    subpool_barcode_column = 1,
    umi_count_column=3,
    read_count_column=4,
  log: '../data/linked_barcode_counts/{illumina_sample_id}.log'
  container: 'docker://ghcr.io/rasilab/python:1.0.0'
  shell:
    """
    awk '
      BEGIN {{
        FS=",";
        OFS=",";
        SUBSEP=",";
        # column to look for linked barcode
        linked_barcode_column={params.offset_to_linked_barcode}+{params.linkage_ref_barcode_number};
        print "subpool_barcode,insert_num,barcode_num,umi_count,read_count"
      }}
      NR == FNR {{
        linkage[$linked_barcode_column] = ${params.linked_insert_num_column}","${params.linked_barcode_num_column};
        next
      }};
      NR > FNR {{
        if (${params.barcode_seq_column} in linkage) {{
          print ${params.subpool_barcode_column},linkage[${params.barcode_seq_column}],${params.umi_count_column},${params.read_count_column}
        }}
      }}
      ' \
      {input.barcode_linkage_file} {input.raw_count_file} 1> {output.linked_count_file} 2> {log}
    """


rule extract_subpool_barcodes_and_count:
  """Extract only reads with right subpool barcodes (upto Hamming distance of 2) and add up UMIs and reads per barcode
  """
  input:
    linked_count_file =  lambda w: f"../data/linked_barcode_counts/{sample_annotations.loc[w.sample_id, 'illumina_sample_id']}.csv",
    notebook = 'get_subpool_counts.ipynb',
  output: '../data/subpool_barcode_counts/{sample_id}.csv'
  params:
    subpool_barcode = lambda w: sample_annotations.loc[w.sample_id, 'subpool_barcode'],
  log:
    '../data/subpool_barcode_counts/{sample_id}.log'
  container: 'docker://ghcr.io/rasilab/r:1.0.0'
  shell:
    """
    jupyter nbconvert --to script --ExecutePreprocessor.kernel_name=ir {input.notebook}
    notebook={input.notebook}
    script="${{notebook/.ipynb/.r}}"
    Rscript ${{script}} {params.subpool_barcode} {wildcards.sample_id} {input.linked_count_file} {output} &> {log}
    """


rule count_reads_umis_barcodes_per_insert:
  """Count reads, UMIs and barcodes per insert

  Include only barcodes with UMI count >= umi_cutoff
  """
  input:
    subpool_count_file = '../data/subpool_barcode_counts/{sample_id}.csv',
  output:
    insert_count_file = '../data/insert_counts/{sample_id}.csv',
  params:
    umi_cutoff = umi_cutoff,
    is_grna = lambda w: int('grna' in sample_annotations.loc[w.sample_id, 'sample_name']),
  log:
    '../data/insert_counts/{sample_id}.log'
  container: 'docker://ghcr.io/rasilab/python:1.0.0'
  shell:
    """
    awk '
    BEGIN {{
      FS=",";
      OFS=",";
      print "insert_num,barcode_count,umi_count,read_count"
    }}
    NR > 1 {{
      if (!{params.is_grna}) {{
        if ($3 >= {params.umi_cutoff}) {{
          barcode_count[$1]++;
          umi_count[$1] += $3;
          read_count[$1] += $4
        }}
      }}
      else {{
        # no umi for gRNA samples, so treat read count as umi count
        if ($4 >= {params.umi_cutoff}) {{
          barcode_count[$1]++;
          umi_count[$1] += $4;
          read_count[$1] += $4
        }}
      }}
    }}
    END {{
      for (e in barcode_count) {{
        print e,barcode_count[e],umi_count[e],read_count[e] | "sort -t, -k3nr"
      }}
    }}
    ' {input.subpool_count_file} 1> {output.insert_count_file} 2> {log}
    """


rule count_reads_umis_barcodes_per_insert_random_partition:
  """Count reads, UMIs and barcodes per insert assigned to two different groups

  Include only barcodes with UMI count >= umi_cutoff
  """
  input:
    subpool_count_file = '../data/subpool_barcode_counts/{sample_id}.csv',
  output:
    insert_count_file = '../data/insert_counts_random_partition/{sample_id}.csv',
  params:
    umi_cutoff = umi_cutoff,
    is_grna = lambda w: int('grna' in sample_annotations.loc[w.sample_id, 'sample_name']),
  log:
    '../data/insert_counts_random_partition/{sample_id}.log'
  container: 'docker://ghcr.io/rasilab/python:1.0.0'
  shell:
    """
    awk '
    BEGIN {{
      FS=",";
      OFS=",";
      SUBSEP=",";
      print "insert_num,barcode_group,barcode_count,umi_count,read_count"
    }}
    NR > 1 {{
      if (!{params.is_grna}) {{
        if ($3 >= {params.umi_cutoff}) {{
          barcode_count[$1,$5]++;
          umi_count[$1,$5] += $3;
          read_count[$1,$5] += $4
        }}
      }}
      else {{
        # no umi for gRNA samples, so treat read count as umi count
        if ($4 >= {params.umi_cutoff}) {{
          barcode_count[$1,$5]++;
          umi_count[$1,$5] += $4;
          read_count[$1,$5] += $4
        }}
      }}
    }}
    END {{
      for (e in barcode_count) {{
        print e,barcode_count[e],umi_count[e],read_count[e] | "sort -t, -k3nr"
      }}
    }}
    ' {input.subpool_count_file} 1> {output.insert_count_file} 2> {log}
    """


rule extract_library_statistics:
  """Extract library statistics from barcode-UMI counts
  """
  input:
    raw_count_file =  lambda w: f"../data/barcode_umi_and_read_counts/{sample_annotations.loc[w.sample_id, 'illumina_sample_id']}.csv",
    subpool_count_file = '../data/subpool_barcode_counts/{sample_id}.csv',
    notebook = 'calculate_library_statistics.ipynb',
  output:
    library_stats = '../data/library_statistics/{sample_id}.csv'
  params:
    umi_cutoff = umi_cutoff
  log:
    '../data/library_statistics/{sample_id}.log'
  container: 'docker://ghcr.io/rasilab/r:1.0.0'
  shell:
    """
    jupyter nbconvert --to script --ExecutePreprocessor.kernel_name=ir {input.notebook}
    notebook={input.notebook}
    script="${{notebook/.ipynb/.r}}"
    Rscript ${{script}} {input.raw_count_file} {input.subpool_count_file} {params.umi_cutoff} {output.library_stats} &> {log}
    """


rule get_mageck_input:
  """Prepare input for MaGeCK
  """
  input:
    treatment_file = lambda w: f"../data/insert_counts/{sample_annotations.index[sample_annotations['sample_name'] == w.treatment].tolist()[0]}.csv",
    control_file = lambda w: f"../data/insert_counts/{sample_annotations.index[sample_annotations['sample_name'] == w.control].tolist()[0]}.csv",
    insert_annotations = "../../rbp_sgrna_barcode_linkage/annotations/insert_annotations.csv",
    notebook = 'get_mageck_input.ipynb',
  output:
    mageck_inputs = '../data/mageck/{treatment}_vs_{control}/input_counts.tsv'
  params:
    umi_cutoff = umi_cutoff,
  container: 'docker://ghcr.io/rasilab/r:1.0.0'
  shell:
    """
    jupyter nbconvert --to script --ExecutePreprocessor.kernel_name=ir {input.notebook}
    notebook={input.notebook}
    script="${{notebook/.ipynb/.r}}"
    Rscript ${{script}} {input.treatment_file} {input.control_file} {input.insert_annotations} {params.umi_cutoff} {output.mageck_inputs} 
    """


rule get_mageck_input_random_partition:
  """Prepare input for MaGeCK
  """
  input:
    treatment_file = lambda w: f"../data/insert_counts_random_partition/{sample_annotations.index[sample_annotations['sample_name'] == w.treatment].tolist()[0]}.csv",
    control_file = lambda w: f"../data/insert_counts_random_partition/{sample_annotations.index[sample_annotations['sample_name'] == w.control].tolist()[0]}.csv",
    insert_annotations = "../../rbp_sgrna_barcode_linkage/annotations/insert_annotations.csv",
    notebook = 'get_mageck_input_random_partition.ipynb',
  output:
    mageck_inputs = '../data/mageck_random_partition/{treatment}_vs_{control}/input_counts.tsv'
  params:
    umi_cutoff = umi_cutoff,
  container: 'docker://ghcr.io/rasilab/r:1.0.0'
  shell:
    """
    jupyter nbconvert --to script --ExecutePreprocessor.kernel_name=ir {input.notebook}
    notebook={input.notebook}
    script="${{notebook/.ipynb/.r}}"
    Rscript ${{script}} {input.treatment_file} {input.control_file} {input.insert_annotations} {params.umi_cutoff} {output.mageck_inputs} 
    """


rule run_mageck:
  """Run MaGeCK to identify genes and sgRNAs with significant fold changes between reporter of interest and control reporter
  """
  input: 
    counts = '../data/mageck/{treatment}_vs_{control}/input_counts.tsv'
  output:
    mageck_gene_summary = '../data/mageck/{treatment}_vs_{control}/mageck.gene_summary.tsv',
    mageck_sgrna_summary = '../data/mageck/{treatment}_vs_{control}/mageck.sgrna_summary.tsv',
  params:
    output_prefix = '../data/mageck/{treatment}_vs_{control}/mageck'
  container: 'docker://ghcr.io/rasilab/mageck:0.5.9'
  log: '../data/mageck/{treatment}_vs_{control}/mageck_output.log'
  shell:
    """
    mageck test -k {input.counts} -t treatment -c control \
    --sort-criteria pos \
    -n {params.output_prefix} \
    --additional-rra-parameters "--min-number-goodsgrna 3" \
    &> {log}
    mv $(echo {output.mageck_gene_summary} | sed 's/.tsv/.txt/') {output.mageck_gene_summary}
    mv $(echo {output.mageck_sgrna_summary} | sed 's/.tsv/.txt/') {output.mageck_sgrna_summary}
    """


rule run_mageck_random_partition:
  """Run MaGeCK to identify genes and sgRNAs with significant fold changes between reporter of interest and control reporter
  """
  input: '../data/mageck_random_partition/{treatment}_vs_{control}/input_counts.tsv'
  output:
    mageck_gene_summary = '../data/mageck_random_partition/{treatment}_vs_{control}/mageck.gene_summary.tsv',
    mageck_sgrna_summary = '../data/mageck_random_partition/{treatment}_vs_{control}/mageck.sgrna_summary.tsv',
  params:
    output_prefix = '../data/mageck_random_partition/{treatment}_vs_{control}/mageck'
  container: 'docker://ghcr.io/rasilab/mageck:0.5.9'
  log: '../data/mageck/{treatment}_vs_{control}/mageck_output.log'
  shell:
    """
    mageck test -k {input} -t treatment -c control \
    --sort-criteria pos \
    -n {params.output_prefix} &> {log}
    mv $(echo {output.mageck_gene_summary} | sed 's/.tsv/.txt/') {output.mageck_gene_summary}
    mv $(echo {output.mageck_sgrna_summary} | sed 's/.tsv/.txt/') {output.mageck_sgrna_summary}
    """


rule combine_mageck_tables:
  """Combine MaGeCK tables for supplementary table generation
  """
  input: 
    mageck_gene_summary = expand('../data/mageck/{sample_name[0]}_vs_{sample_name[1]}/mageck.gene_summary.tsv', 
                           sample_name=list(zip(mageck_comparisons['treatment'], mageck_comparisons['control']))),
    mageck_sgrna_summary = expand('../data/mageck/{sample_name[0]}_vs_{sample_name[1]}/mageck.sgrna_summary.tsv', 
                           sample_name=list(zip(mageck_comparisons['treatment'], mageck_comparisons['control']))),
    notebook = "combine_mageck_data.ipynb"
  output:
    gene_summary = '../data/mageck/gene_summary_table.csv.gz',
    sgrna_summary = '../data/mageck/sgrna_summary_table.csv.gz',
  container: 'docker://ghcr.io/rasilab/r:1.0.0'
  shell:
    """
      jupyter nbconvert --to script --ExecutePreprocessor.kernel_name=ir {input.notebook}
      notebook={input.notebook}
      script="${{notebook/.ipynb/.r}}"
      Rscript ${{script}}
    """