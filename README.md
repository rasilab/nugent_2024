# Decoding post-transcriptional regulatory networks by RNA-linked CRISPR screening in human cells

**Patrick J. Nugent, Heungwon Park, Cynthia L. Wladyka, James N. Yelland, Sayantani Sinha, Katharine Y. Chen, Christine Bynum, Grace Quarterman, Stanley C. Lee, Andrew C. Hsieh & Arvind Rasi Subramaniam**

***Nature Methods*** **22**, 1237–1246 (2025)

**Publication:** [Nature Methods](https://www.nature.com/articles/s41592-025-02702-6)  
**Preprint:** [bioRxiv](https://www.biorxiv.org/content/10.1101/2024.07.25.605204v1.full)  
**Manuscript:** [GitHub-optimized Markdown](manuscript/manuscript.md)  
**Questions?** [GitHub Issues](https://github.com/rasilab/nugent_2024/issues/new/choose)  
**Contact:** [Arvind Rasi Subramaniam](mailto:rasi@fredhutch.org)

**Repository Contents:**
1. [📄 Manuscript](manuscript/)
2. [🔬 Figure Code and Source Data](#-figure-code-and-source-data)
3. [⚙️ How to Reproduce the Analysis](#️-how-to-reproduce-the-analysis)

---

## 🔬 Figure Code and Source Data

This table maps each manuscript figure to its plotting script and source data.

| Figure panel     | Figure | Plotting Code                                                                                          | Source Data                                                                                                        |
| ---------------- | ------ | ------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------ |
| 1b               | [Figure](analysis/flow_cytometry/fig1_eyfp_reporter_sgeyfp/figures/sgyfp_sgfluc_effects_for_validation.pdf) | [Code](analysis/flow_cytometry/fig1_eyfp_reporter_sgeyfp/scripts/plot_fig1_flow.R#L83)                 | [Source data](source_data/figure_1b.csv)                                                                           |
| 1d               | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/compare_mrna_gdna_foldchange_with_time.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_grna_fitness_results.R#L142)                | [Source data](source_data/figure_1d.csv)                                                                           |
| 1e               | [Figure, mRNA](analysis/barcodeseq/rbp_barcode_screens/figures/ntc_total_mrna_foldchange_with_time.pdf) <br> [Figure, gDNA](analysis/barcodeseq/rbp_barcode_screens/figures/gdna_foldchange_with_time.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_grna_fitness_results.R#L169) <br> [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_grna_fitness_results.R#L196) | [Source data, mRNA](source_data/figure_1e_mrna.csv) <br> [Source data, gDNA](source_data/figure_1e_gdna.csv)                   |
| 2a               | [Figure](analysis/polysome_profiling/fig2a_polysome_relic/figures/polysome_profiles.pdf) | [Code](analysis/polysome_profiling/fig2a_polysome_relic/scripts/plot_fig2_polysomes.R#L48)             | [Source data](source_data/figure_2a.csv)                                                                           |
| 2c               | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/polysome_volcano.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.R#L184)                 | [Source data](source_data/figure_2c.csv)                                                                           |
| 2d               | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/polysome_ribosome_groups.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.R#L258)                 | [Source data](source_data/figure_2d.csv)                                                                           |
| 2e               | [Figure, volcano](analysis/barcodeseq/rbp_barcode_screens/figures/polysome_volcano_supernatant.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.R#L220)                 | [Source data, volcano](source_data/figure_2e_volcano.csv) |
| 2f               | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/polysome_translation_groups_2.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.R#L330)                 | [Source data](source_data/figure_2f.csv)                                                                           |
| 2g               | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/polysome_vs_mrna_fitness_all.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.R#L382)                 | [Source data](source_data/figure_2g_s2g_mrna.csv)                                                                  |
| 2h               | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/elongation_vs_initiation_sgrna.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.R#L509)                 | [Source data](source_data/figure_2h.csv)                                                                           |
| 3b               | [Figure](analysis/rnaseq/figures/globin_cvg.png) | [Code](analysis/rnaseq/scripts/make_plasmid_rna_seq_coverage_plots.R#L177)                             | [Source data](source_data/figure_3b.csv)                                                                           |
| 3d               | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/splicing_go_enrichment.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_splicing_results.R#L195)                    | [Source data](source_data/figure_3d.csv)                                                                           |
| 3e               | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/splicing_volcano.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_splicing_results.R#L114)                    | [Source data](source_data/figure_3e.csv)                                                                           |
| 3f               | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/sf3b_lfc.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_splicing_results.R#L151)                    | [Source data](source_data/figure_3f.csv)                                                                           |
| 3g               | [Figure, exon skipping](analysis/rnaseq/figures/skipped_exon_isoform_change.png) <br> [Figure, intron retention](analysis/rnaseq/figures/retained_intron_isoform_change.png) | [Code](analysis/rnaseq/scripts/analyze_exon_skipping.R#L333) <br> [Code](analysis/rnaseq/scripts/analyze_intron_coverage_genome.R#L192) | [Source data, exon skipping](source_data/figure_3g.csv) <br> [Source data, intron retention](source_data/figure_3g.csv) |
| 3h               | [Figure, RPL41](analysis/rnaseq/figures/rpl41_cvg.png) <br> [Figure, RPL24](analysis/rnaseq/figures/rpl24_cvg.png) | [Code](analysis/rnaseq/scripts/make_rna_seq_coverage_plots.R#L257) <br> [Code](analysis/rnaseq/scripts/make_rna_seq_coverage_plots.R#L266) | [Source data, RPL41](source_data/figure_3h.csv) <br> [Source data, RPL24](source_data/figure_3h.csv) |
| 3i               | [Figure, exon skipping](analysis/rnaseq/figures/rpl24_rpl41_exon_skipped_fraction.png) <br> [Figure, intron retention](analysis/rnaseq/figures/rpl24_rpl41_intron_retained_fraction.png) | [Code](analysis/rnaseq/scripts/analyze_exon_skipping.R#L239) <br> [Code](analysis/rnaseq/scripts/analyze_intron_coverage_genome.R#L114) | [Source data, exon skipping](source_data/figure_3i.csv) <br> [Source data, intron retention](source_data/figure_3i.csv) |
| 4b               | [Figure, alphabetic](analysis/barcodeseq/rbp_barcode_screens/figures/nmd_alphabetic.pdf) <br> [Figure, hits jitter](analysis/barcodeseq/rbp_barcode_screens/figures/nmd_hits_jitter.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_nmd_results.R#L121) <br> [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_nmd_results.R#L182) | [Source data](source_data/figure_4b_s4c.csv)                                                                       |
| 4c, 4d           | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/nmd_gcn1_isrib_volcano.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_nmd_results.R#L215)                         | [Source data](source_data/figure_4c_4d.csv)                                                                        |
| 5b               | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/eyfp_harr_volcano.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_eyfp_deopt_harr_results.R#L83)             | [Source data](source_data/figure_5b.csv)                                                                           |
| 5c               | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/gcn1_sgrna_hht.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_eyfp_deopt_harr_results.R#L112)             | [Source data](source_data/figure_5c.csv)                                                                           |
| 5e               | [Figure, scatter](analysis/rnaseq/figures/scatter_plot_hht_treatment_gcn1_vs_fluc.pdf) <br> [Figure, IEG](analysis/rnaseq/figures/ieg_alone_plot_hht_treatment_gcn1_vs_fluc.pdf) | [Code](analysis/rnaseq/scripts/analyze_fold_changes.R#L137) <br> [Code](analysis/rnaseq/scripts/analyze_fold_changes.R#L163) | [Source data, scatter](source_data/figure_5e_scatter.csv) <br> [Source data, IEG](source_data/figure_5e_ieg.csv) |
| 5f               | [Figure](analysis/riboseq/figures/riboseq_metadensity.pdf) | [Code](analysis/riboseq/scripts/analyze_transcriptome_coverage.R#L104)                                   | [Source data](source_data/figure_5f.csv)                                                                           |
| 5g               | [Figure, JUN](analysis/riboseq/figures/jun_riboseq.pdf) <br> [Figure, MYC](analysis/riboseq/figures/myc_riboseq.pdf) | [Code](analysis/riboseq/scripts/analyze_transcriptome_coverage.R#L168) <br> [Code](analysis/riboseq/scripts/analyze_transcriptome_coverage.R#L217) | [Source data, JUN](source_data/figure_5g_jun.csv) <br> [Source data, MYC](source_data/figure_5g_myc.csv) |
| Extended data 1a | [Figure](analysis/flow_cytometry/figs1a_eyfp_reporter_sgeyfp_u2os_293t/figures/sgyfp_sgfluc_effects_for_validation_s1a.pdf) | [Code](analysis/flow_cytometry/figs1a_eyfp_reporter_sgeyfp_u2os_293t/scripts/plot_figs1a_flow.R#L94)   | [Source data](source_data/figure_s1a.csv)                                                                          |
| Extended data 1b | [Figure](analysis/flow_cytometry/figs1b_integration_efficiency_u2os_293t/figures/bfp_vs_rfp_s1b.pdf) | [Code](analysis/flow_cytometry/figs1b_integration_efficiency_u2os_293t/scripts/plot_figs1b_flow.R#L69) | [Source data](source_data/figure_s1b.csv)                                                                          |
| Extended data 1d | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/pdf_sgrna_umi_counts.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_grna_fitness_results.R#L86)                | [Source data](source_data/figure_s1d.csv)                                                                          |
| Extended data 1e | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/cdf_sgrna_n_barcodes.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_grna_fitness_results.R#L68)                | [Source data](source_data/figure_s1e.csv)                                                                          |
| Extended data 1f | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/mageck_fitness_compare_barcode_scatter.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/compare_barcode_partitions.R#L188)               | [Source data](source_data/figure_s1f.csv)                                                                          |
| Extended data 2a | [Figure, histogram](analysis/barcodeseq/rbp_barcode_screens/figures/mageck_polysome_histogram.pdf) <br> [Figure, correlation](analysis/barcodeseq/rbp_barcode_screens/figures/mageck_polysome_correlation.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.R#L148) <br> [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.R#L134) | [Source data](source_data/figure_s2a_histo.csv) <br> [Source data](source_data/figure_s2a_scatter.csv)             |
| Extended data 2b | Figure not generated | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.R#L516)                 | [Source data](source_data/figure_s2b.csv)                                                                          |
| Extended data 2c | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/supernatant_ribosome_groups.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.R#L292)                 | [Source data](source_data/figure_2e_s2c_ribosome_groups.csv)                                                                          |
| Extended data 2d | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/interesting_supernatant_genes.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.R#L112)                 | [Source data](source_data/figure_s2d.csv)                                                                          |
| Extended data 2e | [Figure](analysis/polysome_profiling/figs2e_polysome_relic_hits/figures/normalized_trace_with_PbyM_s2e.pdf) | [Code](analysis/polysome_profiling/figs2e_polysome_relic_hits/scripts/plot_figs2_polysomes.R#L104)      | [Source data](source_data/figure_s2e.csv)                                                                          |
| Extended data 2f | [Figure, gDNA](analysis/barcodeseq/rbp_barcode_screens/figures/polysome_vs_grna_fitness_all.pdf) <br> [Figure, mRNA](analysis/barcodeseq/rbp_barcode_screens/figures/polysome_vs_mrna_fitness_all.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.R#L401) <br> [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.R#L382) | [Source data, gDNA](source_data/figure_s2f_gdna.csv) <br> [Source data, mRNA](source_data/figure_s2f_mrna.csv)                 |
| Extended data 2g | [Figure, mRNA](analysis/barcodeseq/rbp_barcode_screens/figures/polysome_vs_mrna_fitness.pdf) <br> [Figure, gDNA](analysis/barcodeseq/rbp_barcode_screens/figures/polysome_vs_grna_fitness.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.R#L442) <br> [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.R#L483) | [Source data, mRNA](source_data/figure_2g_s2g_mrna.csv) <br> [Source data, gDNA](source_data/figure_s2g_gdna.csv) |
| Extended data 3a | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/splicing_n_hits.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_splicing_results.R#L80)                    | [Source data](source_data/figure_s3a.csv)                                                                          |
| Extended data 3b | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/mageck_splicing_scatter.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/compare_barcode_partitions.R#L91)               | [Source data](source_data/figure_s3b.csv)                                                                          |
| Extended data 3c | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/mageck_splicing_correlation.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/compare_barcode_partitions.R#L127)               | [Source data](source_data/figure_s3c.csv)                                                                          |
| Extended data 3d | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/sf3b_fitness.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_grna_fitness_results.R#L230)                | [Source data](source_data/figure_s3d.csv)                                                                          |
| Extended data 4a | [Figure](analysis/qpcr/figs4_nmd_reporter_validation/figures/mcherry_normalized_ct_values_inverted_s4a.pdf) | [Code](analysis/qpcr/figs4_nmd_reporter_validation/scripts/plot_figs4_qpcr.R#L68)                      | [Source data](source_data/figure_s4a.csv)                                                                          |
| Extended data 4c | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/nmd_volcano.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_nmd_results.R#L134)                         | [Source data](source_data/figure_s4c.csv)                                                                          |
| Extended data 4d | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/ptc_mrna_eif2_eif3_eif4.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_nmd_results.R#L159)                         | [Source data](source_data/figure_s4d.csv)                                                                          |
| Extended data 4e | [Figure](analysis/barcodeseq/rbp_barcode_screens/figures/eif_fitness_mrna.pdf) | [Code](analysis/barcodeseq/rbp_barcode_screens/scripts/plot_grna_fitness_results.R#L278)                | [Source data](source_data/figure_s4e.csv)                                                                          |
| Extended data 5c | [Figure](analysis/qpcr/figs5c_u937_gcn1_hht/figures/egr1_jun_mrna_levels_s5c.pdf) | [Code](analysis/qpcr/figs5c_u937_gcn1_hht/scripts/plot_figs5c_qpcr.R#L80)                              | [Source data](source_data/figure_s5c.csv)                                                                          |
| Extended data 5d | [Figure](analysis/qpcr/figs5d_zaki_gcn1_hht/figures/egr1_jun_mrna_levels_s5d.pdf) | [Code](analysis/qpcr/figs5d_zaki_gcn1_hht/scripts/plot_fig_s5d_qpcr.R#L80)                             | [Source data](source_data/figure_s5d.csv)                                                                          |
| Extended data 5e | [Figure](analysis/polysome_profiling/figs5e_hht_gcn1_mnase/figures/polysome_profiles_mnase.pdf) | [Code](analysis/polysome_profiling/figs5e_hht_gcn1_mnase/scripts/plot_figs5_polysomes.ipynb)           | [Source data](source_data/figure_s5e.csv)                                                                          |

---

## ⚙️ How to Reproduce the Analysis

This repository follows best practices for computational reproducibility in the life sciences, as described in [Grüning *et al.*, 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6263957/).

### 🔧 System Requirements

**Container Runtime (Required):**
Install Apptainer system-wide for container execution:
```bash
# Ubuntu/Debian
sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt update && sudo apt install -y apptainer
```

**Other Dependencies:**
- Git (for cloning repository)
- wget (for downloading conda installer) 
- At least 50GB free disk space

**Alternative Setups:**
- **Fred Hutch users:** Use `module load Singularity snakemake` instead
- **Docker users:** Can substitute Docker for Apptainer (requires modifying submit_local.sh)

### 🚀 Run the Full Analysis

**Local execution:**
```bash
bash run_everything.sh
```

**SLURM cluster execution:**
```bash
SUBMIT_SCRIPT=submit_cluster.sh bash run_everything.sh
```

The script automatically:
* Installs local Conda and Snakemake
* Downloads all required containers from [rasilab GitHub Packages](https://github.com/orgs/rasilab/packages)
* Downloads high-throughput sequencing data from NCBI SRA
* Processes barcode-sgRNA linkage and barcode sequencing data
* Analyzes RNA-seq and ribosome profiling data
* Generates all manuscript figures and source data

### 🛠️ Advanced Usage

**Partial analysis:**
* Copy specific sections from [run_everything.sh](run_everything.sh) for selected analyses

**Individual workflow execution:**
* Local: `cd analysis/<workflow>/scripts && bash submit_local.sh -s run_analysis.smk`
* SLURM cluster: `cd analysis/<workflow>/scripts && bash submit_cluster.sh -s run_analysis.smk`
* Modify [`analysis/cluster.yaml`](./analysis/cluster.yaml) for different SLURM cluster configurations

**Interactive analysis:**
* Use [r_python container](https://github.com/rasilab/r_python/pkgs/container/r_python) with [VSCode](https://rasilab.github.io/docs/software/how_to_create_and_use_containers/)
* See [container documentation](https://rasilab.github.io/docs/software/how_to_create_and_use_containers/) for optimization tips
