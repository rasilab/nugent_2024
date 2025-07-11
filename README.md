# Decoding post-transcriptional regulatory networks by RNA-linked CRISPR screening in human cells


**Publication:** [Nature Methods](https://www.nature.com/articles/s41592-025-02702-6)  
**Preprint:** [bioRxiv](https://www.biorxiv.org/content/10.1101/2024.07.25.605204v1.full)  
**Questions?** [GitHub Issues](https://github.com/rasilab/nugent_2024/issues/new/choose)  
**Contact:** [Arvind Rasi Subramaniam](mailto:rasi@fredhutch.org), [Patrick Nugent](mailto:pnugent@fredhutch.org)

**Repository Contents:**
1. [🔬 Figure Code and Source Data](#-figure-code-and-source-data)
2. [⚙️ How to Reproduce the Analysis](#️-how-to-reproduce-the-analysis)
   1. [✅ Requirements](#-requirements)
   2. [📦 Containers](#-containers)
   3. [🚀 Run the Full Analysis](#-run-the-full-analysis)
   4. [🧪 Interactive Analysis in Jupyter](#-interactive-analysis-in-jupyter)

---

## 🔬 Figure Code and Source Data

This table maps each manuscript figure to its plotting script and source data.

| Figure panel     | Figure | Plotting Code                                                                                          | Source Data                                                                                                        |
| ---------------- | ------ | ------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------ |
| 1b               | [Figure](analysis/flow_cytometry/fig1_eyfp_reporter_sgeyfp/figures/sgyfp_sgfluc_effects_for_validation.pdf) | [Code](analysis/flow_cytometry/fig1_eyfp_reporter_sgeyfp/scripts/plot_fig1_flow.R)                 | [Source data](source_data/figure_1b.csv)                                                                           |
| 1d               | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/compare_mrna_gdna_foldchange_with_time.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_grna_fitness_results.R)                | [Source data](source_data/figure_1d.csv)                                                                           |
| 1e               | [Figure, mRNA](analysis/barcodeseq/rbp_barcode_screens/figures/ntc_total_mrna_foldchange_with_time.pdf) <br> [Figure, gDNA](analysis/barcodeseq/rbp_barcode_screens/figures/gdna_foldchange_with_time.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_grna_fitness_results.R)                | [Source data, mRNA](source_data/figure_1e_mrna.csv) <br> [Source data, gDNA](source_data/figure_1e_gdna.csv)                   |
| 2a               | [Figure](analysis/polysome_profiling/fig2a_polysome_relic/figures/polysome_profiles.pdf) | [Code](analysis/polysome_profiling/fig2a_polysome_relic/scripts/plot_fig2_polysomes.R#L35)             | [Source data](source_data/figure_2a.csv)                                                                           |
| 2c               | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/polysome_volcano.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.R#L184)                 | [Source data](source_data/figure_2c.csv)                                                                           |
| 2d               | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/polysome_ribosome_groups.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.R#L258)                 | [Source data](source_data/figure_2d.csv)                                                                           |
| 2e               | [Figure, volcano](analysis/barcodeseq/rbp_barcode_screens/figures/polysome_volcano_supernatant.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.R#L220)                 | [Source data, volcano](source_data/figure_2e_volcano.csv) |
| 2e               | [Figure, ribosome groups](analysis/barcodeseq/rbp_barcode_screens/figures/supernatant_ribosome_groups.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.R#L292)                 | [Source data, ribosome groups](source_data/figure_2e_s2c_ribosome_groups.csv) |
| 2f               | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/polysome_translation_groups_2.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.R#L330)                 | [Source data](source_data/figure_2f.csv)                                                                           |
| 2g               | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/polysome_vs_mrna_fitness.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.R#L442)                 | [Source data](source_data/figure_2g_s2g_mrna.csv)                                                                  |
| 2h               | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/elongation_vs_initiation_sgrna.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.R#L509)                 | [Source data](source_data/figure_2h.csv)                                                                           |
| 3b               | [Figure](analysis/rnaseq/figures/globin_cvg.png) | [Code](analysis/rnaseq/scripts/make_plasmid_rna_seq_coverage_plots.R#L177)                             | [Source data](source_data/figure_3b.csv)                                                                           |
| 3d               | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/splicing_go_enrichment.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_splicing_results.R#L195)                    | [Source data](source_data/figure_3d.csv)                                                                           |
| 3e               | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/splicing_volcano.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_splicing_results.R#L114)                    | [Source data](source_data/figure_3e.csv)                                                                           |
| 3f               | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/sf3b_lfc.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_splicing_results.R#L151)                    | [Source data](source_data/figure_3f.csv)                                                                           |
| 3g               | [Figure, exon skipping](analysis/rnaseq/figures/skipped_exon_isoform_change.png) <br> [Figure, intron retention](analysis/rnaseq/figures/retained_intron_isoform_change.png) | [Code](analysis/rnaseq/scripts/analyze_exon_skipping.R#L333) <br> [Code](analysis/rnaseq/scripts/analyze_intron_coverage_genome.R#L192) | [Source data, exon skipping](source_data/figure_3g.csv) <br> [Source data, intron retention](source_data/figure_3g.csv) |
| 3h               | [Figure, RPL41](analysis/rnaseq/figures/rpl41_cvg.png) <br> [Figure, RPL24](analysis/rnaseq/figures/rpl24_cvg.png) | [Code](analysis/rnaseq/scripts/make_rna_seq_coverage_plots.R#L257) <br> [Code](analysis/rnaseq/scripts/make_rna_seq_coverage_plots.R#L266) | [Source data, RPL41](source_data/figure_3h.csv) <br> [Source data, RPL24](source_data/figure_3h.csv) |
| 3i               | [Figure, exon skipping](analysis/rnaseq/figures/rpl24_rpl41_exon_skipped_fraction.png) <br> [Figure, intron retention](analysis/rnaseq/figures/rpl24_rpl41_intron_retained_fraction.png) | [Code](analysis/rnaseq/scripts/analyze_exon_skipping.R#L239) <br> [Code](analysis/rnaseq/scripts/analyze_intron_coverage_genome.R#L114) | [Source data, exon skipping](source_data/figure_3i.csv) <br> [Source data, intron retention](source_data/figure_3i.csv) |
| 4b               | [Figure, alphabetic](analysis/barcodeseq/rbp_barcode_screens/figures/nmd_alphabetic.pdf) <br> [Figure, hits jitter](analysis/barcodeseq/rbp_barcode_screens/figures/nmd_hits_jitter.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_nmd_results.R#L121)                         | [Source data](source_data/figure_4b_s4c.csv)                                                                       |
| 4c, 4d           | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/nmd_gcn1_isrib_volcano.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_nmd_results.R#L134)                         | [Source data](source_data/figure_4c_4d.csv)                                                                        |
| 5b               | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/eyfp_harr_volcano.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_eyfp_deopt_harr_results.R#L83)             | [Source data](source_data/figure_5b.csv)                                                                           |
| 5c               | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/gcn1_sgrna_hht.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_eyfp_deopt_harr_results.R#L112)             | [Source data](source_data/figure_5c.csv)                                                                           |
| 5f               | [Figure](analysis/riboseq/figures/riboseq_metadensity.pdf) | [Code](analysis/riboseq/scripts/analyze_transcriptome_coverage.R#L94)                                   | [Source data](source_data/figure_5f.csv)                                                                           |
| 5g               | [Figure, JUN](analysis/riboseq/figures/jun_riboseq.pdf) <br> [Figure, MYC](analysis/riboseq/figures/myc_riboseq.pdf) | [Code](analysis/riboseq/scripts/analyze_transcriptome_coverage.R#L149) <br> [Code](analysis/riboseq/scripts/analyze_transcriptome_coverage.R#L196) | [Source data, JUN](source_data/figure_5g_jun.csv) <br> [Source data, MYC](source_data/figure_5g_myc.csv) |
| Extended data 1a | [Code](analysis/flow_cytometry/figs1a_eyfp_reporter_sgeyfp_u2os_293t/scripts/plot_figs1a_flow.ipynb)   | [Source data](source_data/figure_s1a.csv)                                                                          |
| Extended data 1b | [Code](analysis/flow_cytometry/figs1b_integration_efficiency_u2os_293t/scripts/plot_figs1b_flow.ipynb) | [Source data](source_data/figure_s1b.csv)                                                                          |
| Extended data 1d | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_grna_fitness_results.R)                | [Source data](source_data/figure_s1d.csv)                                                                          |
| Extended data 1e | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_grna_fitness_results.R)                | [Source data](source_data/figure_s1e.csv)                                                                          |
| Extended data 1f | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/compare_barcode_partitions.R)               | [Source data](source_data/figure_s1f.csv)                                                                          |
| Extended data 2a | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.R)                 | [Source data](source_data/figure_s2a_histo.csv) <br> [Source data](source_data/figure_s2a_scatter.csv)             |
| Extended data 2b | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.R)                 | [Source data](source_data/figure_s2b.csv)                                                                          |
| Extended data 2d | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.R)                 | [Source data](source_data/figure_s2d.csv)                                                                          |
| Extended data 2e | [Code](analysis/polysome_profiling/figs2e_polysome_relic_hits/scripts/plot_figs2_polysomes.ipynb)      | [Source data](source_data/figure_s2e.csv)                                                                          |
| Extended data 2f | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.R)                 | [Source data](source_data/figure_s2f_gdna.csv) <br> [Source data](source_data/figure_s2f_mrna.csv)                 |
| Extended data 2g | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.R)                 | [Source data](source_data/figure_s2g_gdna.csv)                                                                     |
| Extended data 3a | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_splicing_results.R)                    | [Source data](source_data/figure_s3a.csv)                                                                          |
| Extended data 3b | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/compare_barcode_partitions.R)               | [Source data](source_data/figure_s3b.csv)                                                                          |
| Extended data 3c | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/compare_barcode_partitions.R)               | [Source data](source_data/figure_s3c.csv)                                                                          |
| Extended data 3d | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_grna_fitness_results.R)                | [Source data](source_data/figure_s3d.csv)                                                                          |
| Extended data 4a | [Code](analysis/qpcr/figs4_nmd_reporter_validation/scripts/plot_figs4_qpcr.ipynb)                      | [Source data](source_data/figure_s4a.csv)                                                                          |
| Extended data 4b | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_nmd_results.R)                         | [Source data](source_data/figure_s4b.csv)                                                                          |
| Extended data 4d | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_nmd_results.R)                         | [Source data](source_data/figure_s4d.csv)                                                                          |
| Extended data 4e | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_grna_fitness_results.R)                | [Source data](source_data/figure_s4e.csv)                                                                          |
| Extended data 5c | [Code](analysis/qpcr/figs5c_u937_gcn1_hht/scripts/plot_figs5c_qpcr.ipynb)                              | [Source data](source_data/figure_s5c.csv)                                                                          |
| Extended data 5d | [Code](analysis/qpcr/figs5d_zaki_gcn1_hht/scripts/plot_fig_s5d_qpcr.ipynb)                             | [Source data](source_data/figure_s5d.csv)                                                                          |
| Extended data 5e | [Code](analysis/polysome_profiling/figs5e_hht_gcn1_mnase/scripts/plot_figs5_polysomes.ipynb)           | [Source data](source_data/figure_s5e.csv)                                                                          |

---

## ⚙️ How to Reproduce the Analysis

This repository follows best practices for computational reproducibility in the life sciences, as described in [Grüning *et al.*, 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6263957/).

### ✅ Requirements

* **Conda (Package Manager)**  
  Install Miniconda (if not already installed):

  ```bash
  # For Linux
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash Miniconda3-latest-Linux-x86_64.sh
  ```

* **Snakemake (Workflow Engine)**  
  Install via Conda:

  ```bash
  conda install -c bioconda snakemake
  ```

* **Singularity (Container Runtime)**  
  Install via Conda:

  ```bash
  conda install -c conda-forge singularity
  ```

Note that Singularity has been renamed to Apptainer and you will see references to both names in online documentation.

Fred Hutch cluster users do not need Conda and can load Snakemake and Singularity as follows. See [here](https://rasilab.github.io/docs/software/how_to_create_and_use_containers/) for more details.

```bash
module load Singularity snakemake
```

---

### 📦 Containers

All tools are packaged as Docker containers:
👉 [rasilab GitHub Packages](https://github.com/orgs/rasilab/packages)

These are pulled automatically when running workflows with Snakemake.

Subramaniam lab users should follow the instructions [here](https://rasilab.github.io/docs/software/how_to_create_and_use_containers/#how-to-use-singularity-containers-in-snakemake-workflows-on-the-fred-hutch-cluster) to symbolically link the `rasilab` container directory rather than downloading the containers each time.

---

### 🚀 Run the Full Analysis

Ensure required paths are mounted using `--bind` in:

* [`analysis/submit_cluster.sh`](./analysis/submit_cluster.sh)
* [`analysis/submit_local.sh`](./analysis/submit_local.sh)

Then run:

```bash
sh run_everything.sh
```

[run_everything.sh](./run_everything.sh) will:

* Download raw high-throughput sequencing data from SRA and convert to FASTQ
* Process linkage and barcode sequencing 
* Analyze RNA-seq
* Analyze Ribo-seq
* Regenerate all manuscript figure panels and source data (excluding schematics and gel images)

> ℹ️ You can also copy specific lines from `run_everything.sh` to reproduce only selected parts.

---

### 🧪 Interactive Analysis in Jupyter

You can explore the data interactively using the Jupyter-ready container:

👉 [r_python container](https://github.com/rasilab/r_python/pkgs/container/r_python)

To use the above container in VSCode, see instructions [here](https://rasilab.github.io/docs/software/how_to_create_and_use_containers/).
