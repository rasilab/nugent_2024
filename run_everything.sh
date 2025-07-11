#!/bin/bash

#==============================================================================
# MASTER SCRIPT: Generate All Manuscript Figures
#==============================================================================
# 
# This script regenerates all figures for the manuscript:
# "Decoding post-transcriptional regulatory networks by RNA-linked CRISPR 
# screening in human cells" (Nature Methods)
#
# USAGE:
#   bash run_everything.sh
#
# REQUIREMENTS:
#   - R with required packages (tidyverse, rasilabRtemplates, etc.)
#   - Singularity container system for specialized analyses
#   - Data files in appropriate directories
#
# OUTPUT:
#   - Figure PDFs in analysis/*/figures/ directories
#   - Source data CSV files in source_data/ directory
#
#==============================================================================

set -e  # Exit on any error

echo "=== Generating All Manuscript Figures ==="

#==============================================================================
# CONTAINER SETUP
#==============================================================================
echo "=== SETTING UP CONTAINERS ==="

# Create local containers directory if it doesn't exist
mkdir -p containers

# Pull DESeq2 container for RNA-seq analysis
echo "Pulling DESeq2 container for RNA-seq analysis..."
if [ ! -f containers/deseq2_1.38.0.simg ]; then
    singularity pull containers/deseq2_1.38.0.simg oras://ghcr.io/rasilab/deseq2:1.38.0
    echo "✓ DESeq2 container downloaded"
else
    echo "✓ DESeq2 container already available"
fi

echo "✓ Container setup completed!"

#==============================================================================
# FIGURE 1: CRISPR Screen Validation and Fitness Analysis
#==============================================================================
echo "=== GENERATING FIGURE 1 ==="

# Figure 1b: Flow cytometry validation of sgRNA effects
echo "Running Figure 1b - Flow cytometry sgRNA validation"
cd analysis/flow_cytometry/fig1_eyfp_reporter_sgeyfp/scripts
Rscript plot_fig1_flow.R
cd ../../../../

# Figure 1d,1e: sgRNA fitness analysis (gDNA vs mRNA, essential gene depletion)
echo "Running Figure 1d,1e - sgRNA fitness and essential gene analysis"
cd analysis/barcodeseq/rbp_barcode_screens/scripts
Rscript plot_grna_fitness_results.R
cd ../../../../

echo "✓ Figure 1 generation completed!"

#==============================================================================
# FIGURE 3: RNA-seq Splicing Analysis
#==============================================================================
echo "=== GENERATING FIGURE 3 ==="

# Run RNA-seq analysis workflows to generate data and figures
echo "Running RNA-seq genome analysis workflow"
cd analysis/rnaseq/scripts
bash ../submit_local.sh run_analysis.smk
cd ../../../

echo "Running RNA-seq plasmid analysis workflow" 
cd analysis/rnaseq/scripts
bash ../submit_local.sh run_analysis_plasmid.smk
cd ../../../

# Figure 5e: RNA-seq fold change analysis using DESeq2 container
echo "Running Figure 5e - RNA-seq fold change analysis (GCN1 vs FLUC with HHT)"
cd analysis/rnaseq
singularity exec ../../containers/deseq2_1.38.0.simg Rscript scripts/analyze_fold_changes.R
cd ../../

# Figure 3d,3e,3f: Splicing screen results
echo "Running Figure 3d,3e,3f - Splicing screen analysis"
cd analysis/barcodeseq/rbp_barcode_screens/scripts
Rscript plot_splicing_results.R
cd ../../../../

echo "✓ Figure 3 generation completed!"

#==============================================================================
# FIGURE 2: Polysome Profiling Analysis  
#==============================================================================
echo "=== GENERATING FIGURE 2 ==="

# Figure 2a: Polysome profiling traces
echo "Running Figure 2a - Polysome profiling traces"
cd analysis/polysome_profiling/fig2a_polysome_relic/scripts
Rscript plot_fig2_polysomes.R
cd ../../../../

# Figure 2c,2d,2e,2f,2g,2h: Polysome ReLiC data
echo "Running Figure 2c-2h - Polysome profiling analysis"
cd analysis/barcodeseq/rbp_barcode_screens/scripts
Rscript plot_polysome_relic_data.R
cd ../../../../

echo "✓ Figure 2 generation completed!"

#==============================================================================
# FIGURE 4: NMD Analysis
#==============================================================================
echo "=== GENERATING FIGURE 4 ==="

# Figure 4b,4c,4d: NMD screen results
echo "Running Figure 4b,4c,4d - NMD analysis"
cd analysis/barcodeseq/rbp_barcode_screens/scripts
Rscript plot_nmd_results.R
cd ../../../../

echo "✓ Figure 4 generation completed!"

#==============================================================================
# FIGURE 5: Translation Initiation Analysis
#==============================================================================
echo "=== GENERATING FIGURE 5 ==="

# Figure 5b,5c: EYFP deoptimization and Harringtonine results
echo "Running Figure 5b,5c - Translation initiation analysis"
cd analysis/barcodeseq/rbp_barcode_screens/scripts
Rscript plot_eyfp_deopt_harr_results.R
cd ../../../../

# Figure 5f,5g: Ribosome profiling analysis
echo "Running Figure 5f,5g - Ribosome profiling analysis"
cd analysis/riboseq/scripts
echo "  Running Snakemake riboseq preprocessing workflow"
bash ../submit_local.sh run_analysis.smk
echo "  Generating riboseq coverage figures and source data"
Rscript analyze_transcriptome_coverage.R
cd ../../../

echo "✓ Figure 5 generation completed!"

#==============================================================================
# EXTENDED DATA FIGURES
#==============================================================================
echo "=== GENERATING EXTENDED DATA FIGURES ==="

# Extended data figures from barcode partition comparisons
echo "Running Extended data figures - Barcode partition analysis"
cd analysis/barcodeseq/rbp_barcode_screens/scripts
Rscript compare_barcode_partitions.R
cd ../../../../

# Extended data 1a: Flow cytometry sgRNA effects in multiple cell lines
echo "Running Extended data 1a - Flow cytometry sgRNA validation"
cd analysis/flow_cytometry/figs1a_eyfp_reporter_sgeyfp_u2os_293t/scripts
Rscript plot_figs1a_flow.R
cd ../../../../

# Extended data 1b: Integration efficiency analysis
echo "Running Extended data 1b - Integration efficiency flow cytometry"
cd analysis/flow_cytometry/figs1b_integration_efficiency_u2os_293t/scripts
Rscript plot_figs1b_flow.R
cd ../../../../

# Extended data 2e: Polysome profiling of selected hits
echo "Running Extended data 2e - Polysome profiling hits"
cd analysis/polysome_profiling/figs2e_polysome_relic_hits/scripts
Rscript plot_figs2_polysomes.R
cd ../../../../

# Extended data 4a: NMD reporter validation by qPCR
echo "Running Extended data 4a - NMD reporter qPCR validation"
cd analysis/qpcr/figs4_nmd_reporter_validation/scripts
Rscript plot_figs4_qpcr.R
cd ../../../../

# Extended data 5c: GCN1 HHT qPCR in U937 cells
echo "Running Extended data 5c - GCN1 HHT qPCR analysis"
cd analysis/qpcr/figs5c_u937_gcn1_hht/scripts
Rscript plot_figs5c_qpcr.R
cd ../../../../

# Extended data 5d: GCN1 HHT with kinase inhibitors qPCR
echo "Running Extended data 5d - GCN1 HHT kinase inhibitor qPCR"
cd analysis/qpcr/figs5d_zaki_gcn1_hht/scripts
Rscript plot_fig_s5d_qpcr.R
cd ../../../../

echo "✓ Extended data generation completed!"

#==============================================================================
# SUMMARY
#==============================================================================
echo ""
echo "=== PIPELINE SUMMARY ==="
echo "✓ Figure generation pipeline completed successfully!"
echo ""
echo "Generated figures can be found in:"
echo "  - analysis/flow_cytometry/fig1_eyfp_reporter_sgeyfp/figures/"
echo "  - analysis/flow_cytometry/figs1a_eyfp_reporter_sgeyfp_u2os_293t/figures/"
echo "  - analysis/flow_cytometry/figs1b_integration_efficiency_u2os_293t/figures/"
echo "  - analysis/barcodeseq/rbp_barcode_screens/figures/"
echo "  - analysis/rnaseq/figures/"
echo "  - analysis/riboseq/figures/"
echo "  - analysis/polysome_profiling/fig2a_polysome_relic/figures/"
echo "  - analysis/polysome_profiling/figs2e_polysome_relic_hits/figures/"
echo "  - analysis/qpcr/figs4_nmd_reporter_validation/figures/"
echo "  - analysis/qpcr/figs5c_u937_gcn1_hht/figures/"
echo "  - analysis/qpcr/figs5d_zaki_gcn1_hht/figures/"
echo ""
echo "Source data files written to:"
echo "  - source_data/"
echo ""
echo "Containers stored in:"
echo "  - containers/"
echo ""
echo "To view the complete figure table, see README.md"
echo "All done! 🎉"