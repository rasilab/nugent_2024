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
#   - Singularity container system
#   - Data files in appropriate directories
#
# OUTPUT:
#   - Figure PDFs in analysis/*/figures/ directories
#   - Source data CSV files in source_data/ directory
#
#==============================================================================

set -e  # Exit on any error

echo "=== Generating All Manuscript Figures ==="

# Container image for all analyses
CONTAINER="docker://ghcr.io/rasilab/r_python:latest"

#==============================================================================
# DATA PROCESSING WORKFLOWS
#==============================================================================
echo "=== RUNNING DATA PROCESSING WORKFLOWS ==="

# Barcode-sgRNA linkage analysis
echo "Running barcode-sgRNA linkage analysis..."
cd analysis/barcodeseq/rbp_sgrna_barcode_linkage/scripts
bash ../submit_local.sh -s run_analysis.smk
cd ../../../../

# Main barcode sequencing analysis  
echo "Running barcode sequencing analysis..."
cd analysis/barcodeseq/rbp_barcode_screens/scripts
bash ../submit_local.sh -s run_analysis.smk
cd ../../../../

# RNA-seq analysis
echo "Running RNA-seq genome analysis..."
cd analysis/rnaseq/scripts
bash ../submit_local.sh -s run_analysis.smk
cd ../../../

echo "Running RNA-seq plasmid analysis..."
cd analysis/rnaseq/scripts
bash ../submit_local.sh -s run_analysis_plasmid.smk
cd ../../../

# Ribosome profiling analysis
echo "Running ribosome profiling analysis..."
cd analysis/riboseq/scripts
bash ../submit_local.sh -s run_analysis.smk
cd ../../../

#==============================================================================
# FIGURE GENERATION
#==============================================================================
echo "=== GENERATING FIGURES ==="

# Figure 1: CRISPR Screen Validation
echo "Generating Figure 1..."
singularity exec $CONTAINER Rscript analysis/flow_cytometry/fig1_eyfp_reporter_sgeyfp/scripts/plot_fig1_flow.R
singularity exec $CONTAINER Rscript analysis/barcodeseq/rbp_barcode_screens/scripts/plot_grna_fitness_results.R

# Figure 2: Polysome Profiling
echo "Generating Figure 2..."
singularity exec $CONTAINER Rscript analysis/polysome_profiling/fig2a_polysome_relic/scripts/plot_fig2_polysomes.R
singularity exec $CONTAINER Rscript analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.R

# Figure 3: RNA-seq Splicing Analysis
echo "Generating Figure 3..."
singularity exec $CONTAINER Rscript analysis/barcodeseq/rbp_barcode_screens/scripts/plot_splicing_results.R

# Figure 4: NMD Analysis
echo "Generating Figure 4..."
singularity exec $CONTAINER Rscript analysis/barcodeseq/rbp_barcode_screens/scripts/plot_nmd_results.R

# Figure 5: Translation Initiation Analysis
echo "Generating Figure 5..."
singularity exec $CONTAINER Rscript analysis/barcodeseq/rbp_barcode_screens/scripts/plot_eyfp_deopt_harr_results.R
singularity exec $CONTAINER Rscript analysis/rnaseq/scripts/analyze_fold_changes.R
singularity exec $CONTAINER Rscript analysis/riboseq/scripts/analyze_transcriptome_coverage.R

# Extended Data Figures
echo "Generating Extended Data figures..."
singularity exec $CONTAINER Rscript analysis/barcodeseq/rbp_barcode_screens/scripts/compare_barcode_partitions.R
singularity exec $CONTAINER Rscript analysis/flow_cytometry/figs1a_eyfp_reporter_sgeyfp_u2os_293t/scripts/plot_figs1a_flow.R
singularity exec $CONTAINER Rscript analysis/flow_cytometry/figs1b_integration_efficiency_u2os_293t/scripts/plot_figs1b_flow.R
singularity exec $CONTAINER Rscript analysis/polysome_profiling/figs2e_polysome_relic_hits/scripts/plot_figs2_polysomes.R
singularity exec $CONTAINER Rscript analysis/qpcr/figs4_nmd_reporter_validation/scripts/plot_figs4_qpcr.R
singularity exec $CONTAINER Rscript analysis/qpcr/figs5c_u937_gcn1_hht/scripts/plot_figs5c_qpcr.R
singularity exec $CONTAINER Rscript analysis/qpcr/figs5d_zaki_gcn1_hht/scripts/plot_fig_s5d_qpcr.R

#==============================================================================
# SUMMARY
#==============================================================================
echo ""
echo "=== ANALYSIS COMPLETE ==="
echo "✓ All workflows executed successfully!"
echo "✓ All figures generated and saved to analysis/*/figures/"
echo "✓ Source data exported to source_data/"
echo ""
echo "See README.md for complete figure-to-code mapping"
echo "All done! 🎉"