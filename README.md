# Decoding RNA Metabolism by RNA-linked CRISPR Screening in Human Cells

**Patrick J. Nugent**<sup>1,2</sup>, **Heungwon Park**<sup>1</sup>, **Cynthia L. Wladyka**<sup>3</sup>, **Andrew C. Hsieh**<sup>3</sup>, **Arvind Rasi Subramaniam**<sup>1,†</sup>

<sup>1</sup> Basic Sciences Division and Computational Biology Section of the Public Health Sciences Division,
Fred Hutchinson Cancer Center, Seattle, WA 98109, USA, <br/>
<sup>2</sup> Molecular and Cellular Biology Graduate Program, University of Washington, Seattle, WA 98195, USA  <br/>
<sup>3</sup> Human Biology Division, Fred Hutchinson Cancer Center, Seattle, WA 98109, USA

<sup>†</sup> Corresponding author: <rasi@fredhutch.org>

## Abstract

RNAs undergo a complex choreography of metabolic processes in human cells that are regulated by thousands of RNA-associated proteins.
While the effects of individual RNA-associated proteins on RNA metabolism have been extensively characterized, the full complement of regulators for most RNA metabolic events remain unknown.
Here we present a massively parallel RNA-linked CRISPR (ReLiC) screening approach to measure the response of diverse RNA metabolic events to knockout of 2,092 human genes encoding all known RNA-associated proteins.
ReLiC screens highlight widespread yet modular interactions between gene networks regulating splicing, translation, and decay of mRNAs.
When combined with biochemical fractionation of polysomes, ReLiC screening reveals striking pathway-specific coupling between growth fitness and mRNA translation.
Perturbing different stages of translation as well as proteasomal function have differential effects on ribosome occupancy, while perturbing mRNA transcription leaves ribosome occupancy largely intact. 
Isoform-selective ReLiC screens capture differential regulation of intron retention and exon skipping by SF3b complex subunits. 
Modifier screens using ReLiC decipher translational regulators upstream of mRNA decay and uncover a role for the ribosome collision sensor GCN1 during treatment with the anti-leukemic drug homoharringtonine.
Our work demonstrates ReLiC as a versatile platform for discovering and dissecting regulatory principles of human RNA metabolism.


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

## Useful Docker containers for interactive analyses using Jupyter Notebooks

- [R](https://github.com/rasilab/r/pkgs/container/r)
- [Python](https://github.com/rasilab/python/pkgs/container/python)
- [R and Python](https://github.com/rasilab/r_python/pkgs/container/r_python)

## Code to generate figure panels from manuscript

| Figure panels                    | Experiment                                                                                       | Script                                                                                                                                                                                           |
| -------------------------------- | ------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| 1b                               | Flow cytometry of cells expressing sgEYFP/FLUC, Cas9, and EYFP reporter                          | [analysis/flow_cytometry/eyfp_reporter_sgeyfp/scripts/plot_fig1_flow.ipynb](analysis/flow_cytometry/eyfp_reporter_sgeyfp/scripts/plot_fig1_flow.ipynb)                                           |
| 1d,e; Extended data 1d,e, 3d, 4e | Barcode seq of ReLiC fitness screen                                                              | [analysis/barcodeseq/rbp_barcode_screens/scripts/plot_grna_fitness_results.ipynb](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_grna_fitness_results.ipynb)                               |
| 2c-h; Extended data 2a-d,f,g     | Barcode seq of polysome ReLiC screen                                                             | [analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.ipynb](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.ipynb)                                 |
| 3d-f; Extended data 3a           | Barcode seq of ReLiC splicing screen                                                             | [analysis/barcodeseq/rbp_barcode_screens/scripts/plot_splicing_results.ipynb](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_splicing_results.ipynb)                                       |
| 4b-d; Extended data 4b-d         | Barcode seq of ReLiC NMD screen                                                                  | [analysis/barcodeseq/rbp_barcode_screens/scripts/plot_nmd_results.ipynb](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_nmd_results.ipynb)                                                 |
| 5b,c                             | Barcode seq of ReLiC harringtonine screen                                                        | [analysis/barcodeseq/rbp_barcode_screens/scripts/plot_eyfp_deopt_harr_results.ipynb](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_eyfp_deopt_harr_results.ipynb)                         |
| Extended data 1a                 | Flow cytometry of U2OS and 293T expressing sgEYFP/FLUC, Cas9, and EYFP reporter                  | [analysis/flow_cytometry/figs1a_eyfp_reporter_sgeyfp_u2os_293t/scripts/plot_figs1a_flow.ipynb](analysis/flow_cytometry/figs1a_eyfp_reporter_sgeyfp_u2os_293t/scripts/plot_figs1a_flow.ipynb)     |
| Extended data 1b                 | Flow cytometry to measure landing pad integration in U2OS and 293T                               | [analysis/flow_cytometry/figs1b_integration_efficiency_u2os_293t/scripts/plot_figs1b_flow.ipynb](analysis/flow_cytometry/figs1b_integration_efficiency_u2os_293t/scripts/plot_figs1b_flow.ipynb) |
| Extended data 2e                 | Polysome profiles from polysome ReLiC screen hits                                                | [analysis/polysome_profiling/figs2e_polysome_relic_hits/scripts/plot_figs2_polysomes.ipynb](analysis/polysome_profiling/figs2e_polysome_relic_hits/scripts/plot_figs2_polysomes.ipynb)           |
| Extended data 4a                 | RT-qPCR comparison of β-globin PTC/NTC reporter levels                                           | [analysis/qpcr/nmd_reporter_validation/scripts/plot_figs4_qpcr.ipynb](analysis/qpcr/nmd_reporter_validation/scripts/plot_figs4_qpcr.ipynb)                                                       |
| Extended data 5c                 | RT-qPCR measurement of GCN1-dependent EGR1, JUN upregulation during HHT treatment                | [analysis/qpcr/figs5c_u937_gcn1_hht/scripts/plot_figs5c_qpcr.ipynb](analysis/qpcr/figs5c_u937_gcn1_hht/scripts/plot_figs5c_qpcr.ipynb)                                                           |
| Extended data 5d                 | RT-qPCR measurement of GCN1-/HHT-dependent EGR1, JUN upregulation during GCN2 and ZAK inhibition | [analysis/qpcr/figs5d_zaki_gcn1_hht/scripts/plot_fig_s5d_qpcr.ipynb](analysis/qpcr/figs5d_zaki_gcn1_hht/scripts/plot_fig_s5d_qpcr.ipynb)                                                         |
| Extended data 5e                 | Polysome profiles from GCN1-depleted and control cells during HHT treatment                      | [analysis/polysome_profiling/hht_gcn1_mnase/plot_figs4_polysomes.ipynb](analysis/polysome_profiling/hht_gcn1_mnase/plot_figs4_polysomes.ipynb)                                                   |

<!-- To do: 
- update polysome Relic barcode seq analysis to reflect current Fig 2c 
- add analyses for 2a,b,i, 3b,g-i, 5e-g extended data 1f, 3b,c, 5a,f and add to table
- MAKE RUN_EVERYTHING!
-->