# Decoding post-transcriptional regulatory networks by RNA-linked CRISPR screening in human cells

Paper: <https://www.nature.com/articles/s41592-025-02702-6>

Preprint: <https://www.biorxiv.org/content/10.1101/2024.07.25.605204v1.full>

Questions? Email <rasi@fredhutch.org>

</small>

1. [Link to code and source data for manuscript figures](#link-to-code-and-source-data-for-manuscript-figures)
2. [Instructions for running the code](#instructions-for-running-the-code)
3. [Docker container for interactive analyses within Jupyter Notebooks](#docker-container-for-interactive-analyses-within-jupyter-notebooks)


## Link to code and source data for manuscript figures

| Figure panel     | Plotting Code                                                                                          | Source Data                                                                                                        |
| ---------------- | ------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------ |
| 1b               | [Code](analysis/flow_cytometry/fig1_eyfp_reporter_sgeyfp/scripts/plot_fig1_flow.ipynb)                 | [Source data](source_data/figure_1b.csv)                                                                           |
| 1d               | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_grna_fitness_results.ipynb)                | [Source data](source_data/figure_1d.csv)                                                                           |
| 1e               | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_grna_fitness_results.ipynb)                | [Source data](source_data/figure_1e_mrna.csv) <br> [Source data](source_data/figure_1e_gdna.csv)                   |
| Extended data 1d | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_grna_fitness_results.ipynb)                | [Source data](source_data/figure_s1d.csv)                                                                          |
| Extended data 1e | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_grna_fitness_results.ipynb)                | [Source data](source_data/figure_s1e.csv)                                                                          |
| Extended data 3d | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_grna_fitness_results.ipynb)                | [Source data](source_data/figure_s3d.csv)                                                                          |
| Extended data 4e | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_grna_fitness_results.ipynb)                | [Source data](source_data/figure_s4e.csv)                                                                          |
| 2a               | [Code](analysis/polysome_profiling/fig2a_polysome_relic/scripts/plot_fig2_polysomes.ipynb)             | [Source data](source_data/figure_2a.csv)                                                                           |
| 2c               | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.ipynb)                 | [Source data](source_data/figure_2c.csv)                                                                           |
| 2d               | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.ipynb)                 | [Source data](source_data/figure_2d.csv)                                                                           |
| 2e               | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.ipynb)                 | [Source data](source_data/figure_2e_s2c_ribosome_groups.csv) <br> [Source data](source_data/figure_2e_volcano.csv) |
| 2f               | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.ipynb)                 | [Source data](source_data/figure_2f.csv)                                                                           |
| 2g               | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.ipynb)                 | [Source data](source_data/figure_2g_s2g_mrna.csv)                                                                  |
| 2h               | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.ipynb)                 | [Source data](source_data/figure_2h.csv)                                                                           |
| Extended data 2a | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.ipynb)                 | [Source data](source_data/figure_s2a_histo.csv) <br> [Source data](source_data/figure_s2a_scatter.csv)             |
| Extended data 2b | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.ipynb)                 | [Source data](source_data/figure_s2b.csv)                                                                          |
| Extended data 2d | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.ipynb)                 | [Source data](source_data/figure_s2d.csv)                                                                          |
| Extended data 2f | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.ipynb)                 | [Source data](source_data/figure_s2f_gdna.csv) <br> [Source data](source_data/figure_s2f_mrna.csv)                 |
| Extended data 2g | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.ipynb)                 | [Source data](source_data/figure_s2g_gdna.csv)                                                                     |
| 3d               | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_splicing_results.ipynb)                    | [Source data](source_data/figure_3d.csv)                                                                           |
| 3e               | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_splicing_results.ipynb)                    | [Source data](source_data/figure_3e.csv)                                                                           |
| 3f               | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_splicing_results.ipynb)                    | [Source data](source_data/figure_3f.csv)                                                                           |
| Extended data 3a | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_splicing_results.ipynb)                    | [Source data](source_data/figure_s3a.csv)                                                                          |
| 4b               | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_nmd_results.ipynb)                         | [Source data](source_data/figure_4b_s4c.csv)                                                                       |
| 4c               | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_nmd_results.ipynb)                         | [Source data](source_data/figure_4c_4d.csv)                                                                        |
| 4d               | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_nmd_results.ipynb)                         | [Source data](source_data/figure_4c_4d.csv)                                                                        |
| Extended data 4b | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_nmd_results.ipynb)                         | [Source data](source_data/figure_s4b.csv)                                                                          |
| Extended data 4d | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_nmd_results.ipynb)                         | [Source data](source_data/figure_s4d.csv)                                                                          |
| 5b               | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_eyfp_deopt_harr_results.ipynb)             | [Source data](source_data/figure_5b.csv)                                                                           |
| 5c               | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_eyfp_deopt_harr_results.ipynb)             | [Source data](source_data/figure_5c.csv)                                                                           |
| Extended data 1a | [Code](analysis/flow_cytometry/figs1a_eyfp_reporter_sgeyfp_u2os_293t/scripts/plot_figs1a_flow.ipynb)   | [Source data](source_data/figure_s1a.csv)                                                                          |
| Extended data 1b | [Code](analysis/flow_cytometry/figs1b_integration_efficiency_u2os_293t/scripts/plot_figs1b_flow.ipynb) | [Source data](source_data/figure_s1b.csv)                                                                          |
| Extended data 1f | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/compare_barcode_partitions.ipynb)               | [Source data](source_data/figure_s1f.csv)                                                                          |
| Extended data 3b | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/compare_barcode_partitions.ipynb)               | [Source data](source_data/figure_s3b.csv)                                                                          |
| Extended data 3c | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/compare_barcode_partitions.ipynb)               | [Source data](source_data/figure_s3c.csv)                                                                          |
| Extended data 2e | [Code](analysis/polysome_profiling/figs2e_polysome_relic_hits/scripts/plot_figs2_polysomes.ipynb)      | [Source data](source_data/figure_s2e.csv)                                                                          |
| Extended data 4a | [Code](analysis/qpcr/figs4_nmd_reporter_validation/scripts/plot_figs4_qpcr.ipynb)                      | [Source data](source_data/figure_s4a.csv)                                                                          |
| Extended data 5c | [Code](analysis/qpcr/figs5c_u937_gcn1_hht/scripts/plot_figs5c_qpcr.ipynb)                              | [Source data](source_data/figure_s5c.csv)                                                                          |
| Extended data 5d | [Code](analysis/qpcr/figs5d_zaki_gcn1_hht/scripts/plot_fig_s5d_qpcr.ipynb)                             | [Source data](source_data/figure_s5d.csv)                                                                          |
| Extended data 5e | [Code](analysis/polysome_profiling/figs5e_hht_gcn1_mnase/scripts/plot_figs5_polysomes.ipynb)           | [Source data](source_data/figure_s5e.csv)                                                                          |

## Instructions for running the code

- All software necessary for running the code are available as Docker images at https://github.com/orgs/rasilab/packages. These images can be downloaded to your local computer using [Docker](https://www.docker.com/) or [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html). Most distributed clusters will already have Singularity available. You can also download and install [Singularity](https://anaconda.org/conda-forge/singularity) on your local computer using [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).
 
- We use [Snakemake](https://anaconda.org/bioconda/snakemake-minimal) for workflow management. This can be installed using Conda or might be already available in your distributed cluster.

- To reproduce the analysis on a cluster, load Singularity and Snakemake (or activate the Conda environment with these software). Ensure that all necessary folders are mounted using `--bind` in [analysis/submit_cluster.sh](./analysis/submit_cluster.sh) and [analysis/submit_local.sh](./analysis/submit_local.sh). These folder locations will be specific to your computing environment. If the correct location is not mounted, you will get `path not found` error in Snakemake workflows that use Singularity containers.

```
module load Singularity snakemake # for fred hutch cluster
sh run_everything.sh
```

- Typically, you will want to run only parts of the [run_everything.sh](./run_everything.sh) script depending on which analysis or figure you are trying to reproduce. You can paste the corresponding lines from the script onto the command line.

- The [run_everything.sh](./run_everthing.sh) script will:
  - Run linkage sequencing and barcode sequencing [analysis](analysis/barcodeseq)
  - Run RNA seq [analysis](analysis/rnaseq)
  - Run Ribo seq [analysis](analysis/riboseq)
  - Run [code](analysis/run_all_ipynb_scripts.smk) to regenerate figure panels, which will also run the analysis of flow cytometry and polysome profiling data
  - Note that the high throughput sequencing analysis will also automatically download the SRA files and convert them to FASTQ.

## Docker container for interactive analyses within Jupyter Notebooks

- [R and Python](https://github.com/rasilab/r_python/pkgs/container/r_python)
