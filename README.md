# Massively parallel identification of sequence motifs triggering ribosome-associated mRNA quality control

**Katharine Y. Chen**<sup>1,2</sup>, **Heungwon Park**<sup>1</sup>, **Arvind Rasi Subramaniam**<sup>1,†</sup>

<sup>1</sup> Basic Sciences Division and Computational Biology Section of the Public
Health Sciences Division, Fred Hutchinson Cancer Center, Seattle, WA
98109, USA <br/>
<sup>2</sup> Molecular and Cellular Biology Program, University of Washington,
Seattle, WA 98195, USA <br/>

<sup>†</sup> Corresponding author: <rasi@fredhutch.org>

Nucleic Acids Research [10.1093/nar/gkae285](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkae285/7655782)

## Abstract

Decay of mRNAs can be triggered by ribosome slowdown at stretches of rare codons or positively charged amino acids.
However, the full diversity of sequence motifs that trigger co-translational mRNA decay is poorly understood.
To comprehensively identify sequence motifs that trigger mRNA decay, we use a massively parallel reporter assay to measure the effect of all possible combinations of codon pairs on mRNA levels in S. cerevisiae.
In addition to known mRNA-destabilizing sequences, we identify several dipeptide repeats whose translation reduces mRNA levels. 
These include combinations of positively charged and bulky residues, as well as proline-glycine and proline-aspartic acid dipeptide repeats.
Genetic deletion of the ribosome collision sensor Hel2 rescues the mRNA effects of these motifs, suggesting that they trigger ribosome slowdown and activate the ribosome-associated quality control (RQC) pathway.
Deep mutational scanning of an mRNA-destabilizing dipeptide repeat reveals a complex relationship between the charge, bulkiness, and location of amino acid residues in conferring mRNA instability.
Finally, we show that the mRNA effects of codon pairs are predictive of the effects of endogenous sequences.
Our work highlights the complexity of sequence motifs driving co-translational mRNA decay in eukaryotes, and presents a high-throughput approach to dissect their requirements at the codon level.

- [Abstract](#abstract)
- [Instructions for running the code repo](#running-the-code)
- [Data](data/)
  - Includes flow cytometry raw data that is also available at http://flowrepository.org/id/FR-FCM-Z6QH
- [Code for designing endogenous fragments library](analysis/library_design/endogenous_fragments/)
- [Code for processing flow cytometry data and regenerating figure panels](analysis/flow_cytometry/wt_hel2_8xdicodon/scripts)
- [Code for linking barcodes and codon pair inserts in original 8× dicodon library](analysis/barcodeseq/8xdicodon_linkage/scripts/)
- [Code for linking barcodes and codon pair inserts in frameshifted 8× dicodon library](analysis/barcodeseq/frameshifted_8xdicodon_linkage/)
- [Code for linking barcodes and codon pair inserts in small-scale validation library](analysis/barcodeseq/mini_8xdicodon_linkage/scripts/)
- Code for counting barcodes and regenerating figures:
  - [codon pair library wild-type cells](analysis/barcodeseq/wt_mrna_grna/scripts/)
  - [codon pair library *hel2Δ* and *syh1Δ* cells](analysis/barcodeseq/hel2_syh1_mrna_grna/scripts/)
  - [frameshifted codon pair library](analysis/barcodeseq/wt_frameshifted_mrna_grna/scripts/)
  - [small-scale validation library](analysis/barcodeseq/wt_hel2_mini_pool/scripts/)
  - [codon pair library under glucose depletion](analysis/barcodeseq/wt_hel2_no_glucose_mrna_grna/scripts/)
  - [small codon pair library in RQC deletion strains](analysis/barcodeseq/small_8xdicodon_rqcdel_mrna_grna/scripts)
- [Code for counting insert-UMI pairs in (FK)~8~ DMS library in wild-type and *hel2∆* cells](analysis/barcodeseq/wt_hel2_fk8_dms/scripts/)
- [Code for counting insert-UMI pairs in (FK)~8~ DMS library in *upf1∆*_cells](analysis/barcodeseq/upf1_fk8_dms/scripts/)
- [Code for linking barcodes and codon pair inserts in small 8× dicodon library in RQC deletion strains](analysis/barcodeseq/small_8xdicodon_rqcdel_linkage/scripts)
- [Code to regenerate all figure panels](analysis/run_all_ipynb_scripts.smk)

## Running the code
- To run this on a cluster with singularity containers, do:
```
module load singularity # for fred hutch cluster
conda activate snakemake # this is a minimal conda env that has snakemake-minimal and pandas for invoking snakefile
sh run_everything.sh
```

- The ```run_everything.sh``` file will:
  - Download FASTQ files from SRA
  - Run all linkage sequencing, barcode sequencing, and insert sequencing [code](analysis/barcodeseq)
  - Run [code](analysis/library_design/endogenous_fragments/scripts/run_analysis.smk) to design of the endogenous fragments library
  - Run [code](analysis/run_all_ipynb_scripts.smk) to regenerate figure panels
    - This will also run the code to process flow cytometry data

## Docker containers
- [R](https://github.com/rasilab/r/pkgs/container/r)
- [python](https://github.com/rasilab/python/pkgs/container/python)
- [R and python](https://github.com/rasilab/r_python/pkgs/container/r_python)


## Code to generate figure panels from manuscript

| Figure panels        | Experiment                                                                                           | Script                                                                                    |
| -------------------- | ---------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------- |
| 1B, 1C, 1D, 1E, 1F   | Barcode seq of wild-type cells with codon pair library                                               | [analysis/barcodeseq/wt_mrna_grna/scripts/plot_aggregate_effects.ipynb](analysis/barcodeseq/wt_mrna_grna/scripts/plot_aggregate_effects.ipynb)                  |
| 2A, 2B, S1D          | Barcode seq of wild-type cells with codon pair library                                               | [analysis/barcodeseq/wt_mrna_grna/scripts/plot_dipeptide_effects.ipynb](analysis/barcodeseq/wt_mrna_grna/scripts/plot_dipeptide_effects.ipynb)                   |
| 2C, 2D               | Flow cytometry of wild-type cells with individual codon pair inserts                                 | [analysis/flow_cytometry/scripts/plot_figure2_flow.ipynb](analysis/flow_cytometry/scripts/plot_figure2_flow.ipynb)                                 |
| 3A                   | Barcode seq of wild-type cells with codon pair library computationally frameshifted                  | [analysis/barcodeseq/wt_hel2_no_glucose_mrna_grna/scripts/plot_translation_effects.ipynb](analysis/barcodeseq/wt_hel2_no_glucose_mrna_grna/scripts/plot_translation_effects.ipynb) |
| 3B, 3C, S4B          | Barcode seq of wild-type cells with codon pair library during glucose depletion                      | [analysis/barcodeseq/wt_hel2_no_glucose_mrna_grna/scripts/plot_translation_effects.ipynb](analysis/barcodeseq/wt_hel2_no_glucose_mrna_grna/scripts/plot_translation_effects.ipynb) |
| 3E, 3F               | Barcode seq of wild-type cells with -1 frameshifted codon pair library                               | [analysis/barcodeseq/wt_hel2_no_glucose_mrna_grna/scripts/plot_translation_effects.ipynb](analysis/barcodeseq/wt_hel2_no_glucose_mrna_grna/scripts/plot_translation_effects.ipynb) |
| 4C, 4D, 4E, S3C, S6A | Barcode seq of *hel2∆* and *syh1∆* cells with codon pair library                                     | [analysis/barcodeseq/hel2_syh1_mrna_grna/scripts/plot_hel2_syh1_dipeptide_effects.ipynb](analysis/barcodeseq/hel2_syh1_mrna_grna/scripts/plot_hel2_syh1_dipeptide_effects.ipynb)  |
| 5B, 5C, S5C          | Deep mutational scan of (FK)~8~ in wild-type and *hel2∆* cells                                       | [analysis/barcodeseq/wt_hel2_fk8_dms/scripts/plot_variant_effects.ipynb](analysis/barcodeseq/wt_hel2_fk8_dms/scripts/plot_variant_effects.ipynb)                  |
| 6B, 6C, 6D           | Barcode seq of wild-type and *hel2∆* cells with endogenous fragments library                         | [analysis/barcodeseq/endo_frag_mrna_grna/scripts/plot_endogenous_frags.ipynb](analysis/barcodeseq/endo_frag_mrna_grna/scripts/plot_endogenous_frags.ipynb)             |
| S1A                  | Barcode seq of wild-type cells with codon pair library                                               | [analysis/barcodeseq/wt_mrna_grna/scripts/plot_supp_alignment_stats.ipynb](analysis/barcodeseq/wt_mrna_grna/scripts/plot_supp_alignment_stats.ipynb)                |
| S1B, S1C, S1E        | Barcode seq of wild-type cells with codon pair library                                               | [analysis/barcodeseq/wt_mrna_grna/scripts/plot_supplemental_missing_data.ipynb](analysis/barcodeseq/wt_mrna_grna/scripts/plot_supplemental_missing_data.ipynb)           |
| S2A, S2B             | Barcode seq of wild-type with codon pair or mini-pool libraries                                      | [analysis/barcodeseq/wt_mrna_grna/scripts/plot_supp_alignment_stats.ipynb](analysis/barcodeseq/wt_mrna_grna/scripts/plot_supp_alignment_stats.ipynb)                |
| S2C, S2D, S2E        | Flow cytometry of wild-type and *hel2∆* cells with individual codon pair inserts                     | [analysis/barcodeseq/hel2_syh1_mrna_grna/scripts/plot_supp_aln_qc.ipynb](analysis/barcodeseq/hel2_syh1_mrna_grna/scripts/plot_supp_aln_qc.ipynb)                  |
| S3A, S3B             | Summary of wild-type cells with codon pair library during all translation-related conditions         | [analysis/barcodeseq/wt_hel2_no_glucose_mrna_grna/scripts/plot_translation_effects.ipynb](analysis/barcodeseq/wt_hel2_no_glucose_mrna_grna/scripts/plot_translation_effects.ipynb) |
| S4A                  | Wild-type, *hel2∆*, *syh1∆*, *hel2∆/syh1∆*, *cue2∆*, and *xrn1∆* cells with revision codon pair pool | [analysis/barcodeseq/small_8xdicodon_rqcdel_mrna_grna/scripts/plot_dicodon_effects.ipynb](analysis/barcodeseq/small_8xdicodon_rqcdel_mrna_grna/scripts/plot_dicodon_effects.ipynb) |
| S5A, S5B             | Wild-type, *hel2∆*, and *upf1∆* cells with (FK)~8~ DMS library                                       | [analysis/barcodeseq/upf1_fk8_dms/scripts/plot_variant_effects_wt_hel2_upf1_reps.ipynb](analysis/barcodeseq/upf1_fk8_dms/scripts/plot_variant_effects_wt_hel2_upf1_reps.ipynb)   |

