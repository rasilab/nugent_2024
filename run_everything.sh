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

# Figure 3d,3e,3f: Splicing screen results
echo "Running Figure 3d,3e,3f - Splicing screen analysis"
cd analysis/barcodeseq/rbp_barcode_screens/scripts
Rscript plot_splicing_results.R
cd ../../../../

echo "✓ Figure 3 generation completed!"

#==============================================================================
# FUTURE FIGURES (To be added)
#==============================================================================
echo "⚠ Additional figures (2,4,5) will be added in future updates"

# TODO: Add Figure 2 - Polysome profiling analysis
# TODO: Add Figure 4 - NMD analysis
# TODO: Add Figure 5 - Translation initiation analysis

#==============================================================================
# SUMMARY
#==============================================================================
echo ""
echo "=== PIPELINE SUMMARY ==="
echo "✓ Figure generation pipeline completed successfully!"
echo ""
echo "Generated figures can be found in:"
echo "  - analysis/flow_cytometry/fig1_eyfp_reporter_sgeyfp/figures/"
echo "  - analysis/barcodeseq/rbp_barcode_screens/figures/"
echo "  - analysis/rnaseq/figures/"
echo ""
echo "Source data files written to:"
echo "  - source_data/"
echo ""
echo "To view the complete figure table, see README.md"
echo "All done! 🎉"