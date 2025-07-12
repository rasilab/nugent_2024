#!/bin/bash

# Generate all manuscript figures for:
# "Decoding post-transcriptional regulatory networks by RNA-linked CRISPR 
# screening in human cells" (Nature Methods)

set -e  # Exit on any error

echo "=== Generating All Manuscript Figures ==="

# Set up local directories
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOCAL_ENV_DIR="$PROJECT_DIR/.env"
CONDA_DIR="$LOCAL_ENV_DIR/miniconda3"
SINGULARITY_CACHE="$LOCAL_ENV_DIR/singularity_cache"

# Create directories if they don't exist
mkdir -p "$LOCAL_ENV_DIR"
mkdir -p "$SINGULARITY_CACHE"

# Set singularity cache directory
export SINGULARITY_CACHEDIR="$SINGULARITY_CACHE"

# Install conda, snakemake, and singularity if needed

# Check if local conda is installed, install if needed
if [ ! -f "$CONDA_DIR/bin/conda" ]; then
    echo "Installing Miniconda to project directory..."
    CONDA_INSTALLER="$LOCAL_ENV_DIR/Miniconda3-latest-Linux-x86_64.sh"
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O "$CONDA_INSTALLER"
    bash "$CONDA_INSTALLER" -b -p "$CONDA_DIR"
    rm "$CONDA_INSTALLER"
    echo "✓ Miniconda installed to $CONDA_DIR"
else
    echo "✓ Local Conda already available at $CONDA_DIR"
fi

# Setup conda environment
export PATH="$CONDA_DIR/bin:$PATH"
source "$CONDA_DIR/etc/profile.d/conda.sh"

# Create/activate nugent_2024 environment
if ! conda env list | grep -q "^nugent_2024 "; then
    echo "Creating nugent_2024 environment..."
    conda create -n nugent_2024 -c bioconda -c conda-forge snakemake singularity -y
else
    # Install missing packages if needed
    conda activate nugent_2024
    for pkg in "snakemake:bioconda" "singularity:conda-forge"; do
        name="${pkg%:*}"
        channel="${pkg#*:}"
        if ! conda list "$name" &> /dev/null; then
            echo "Installing $name..."
            conda install -c "$channel" "$name" -y
        fi
    done
fi
conda activate nugent_2024
echo "✓ Environment ready"

# Download all required containers

# Define containers to pull
declare -A CONTAINERS=(
    ["r_python_1.3.0"]="docker://ghcr.io/rasilab/r_python:1.3.0"
    ["deseq2_1.38.0"]="docker://ghcr.io/rasilab/deseq2:1.38.0"
    ["sratools_3.0.8"]="docker://ghcr.io/rasilab/sratools:3.0.8"
    ["parallel_fastq_dump_0.6.7"]="docker://ghcr.io/rasilab/parallel_fastq_dump:0.6.7"
    ["python_1.0.0"]="docker://ghcr.io/rasilab/python:1.0.0"
    ["r_1.0.0"]="docker://ghcr.io/rasilab/r:1.0.0"
    ["bowtie_1.3.1"]="docker://ghcr.io/rasilab/bowtie:1.3.1"
    ["bowtie2_2.4.5"]="docker://ghcr.io/rasilab/bowtie2:2.4.5"
    ["cutadapt_4.4"]="docker://ghcr.io/rasilab/cutadapt:4.4"
    ["entrez-direct_16.2"]="docker://ghcr.io/rasilab/entrez-direct:16.2"
    ["mageck_0.5.9"]="docker://ghcr.io/rasilab/mageck:0.5.9"
    ["samtools_1.16.1"]="docker://ghcr.io/rasilab/samtools:1.16.1"
    ["star_2.7.11a"]="docker://ghcr.io/rasilab/star:2.7.11a"
)

# Pull all containers
# Pull containers
for CONTAINER_NAME in "${!CONTAINERS[@]}"; do
    SIF_PATH="$SINGULARITY_CACHE/${CONTAINER_NAME}.sif"
    if [ ! -f "$SIF_PATH" ]; then
        echo "Pulling $CONTAINER_NAME..."
        singularity pull "$SIF_PATH" "${CONTAINERS[$CONTAINER_NAME]}"
    fi
done
echo "✓ All containers ready"

echo "
=== RUNNING DATA PROCESSING WORKFLOWS ==="

# Run Snakemake workflows
run_workflow() {
    local name="$1"
    local path="$2"
    local smk_file="${3:-run_analysis.smk}"
    
    echo "Running $name..."
    (cd "$PROJECT_DIR/$path" && bash submit_local.sh -s "$smk_file")
}

run_workflow "barcode-sgRNA linkage analysis" "analysis/barcodeseq/rbp_sgrna_barcode_linkage/scripts"
run_workflow "barcode sequencing analysis" "analysis/barcodeseq/rbp_barcode_screens/scripts"
run_workflow "RNA-seq genome analysis" "analysis/rnaseq/scripts"
run_workflow "RNA-seq plasmid analysis" "analysis/rnaseq/scripts" "run_analysis_plasmid.smk"
run_workflow "ribosome profiling analysis" "analysis/riboseq/scripts"

echo "
=== GENERATING FIGURES ==="

# Helper function to run R scripts
run_r() {
    local container="${2:-r_python_1.3.0}"
    singularity exec "$SINGULARITY_CACHE/${container}.sif" Rscript "$1"
}

# Main figures
echo "Generating Figure 1 (CRISPR Screen Validation)..."
run_r analysis/flow_cytometry/fig1_eyfp_reporter_sgeyfp/scripts/plot_fig1_flow.R
run_r analysis/barcodeseq/rbp_barcode_screens/scripts/plot_grna_fitness_results.R

echo "Generating Figure 2 (Polysome Profiling)..."
run_r analysis/polysome_profiling/fig2a_polysome_relic/scripts/plot_fig2_polysomes.R
run_r analysis/barcodeseq/rbp_barcode_screens/scripts/plot_polysome_relic_data.R

echo "Generating Figure 3 (RNA-seq Splicing)..."
run_r analysis/barcodeseq/rbp_barcode_screens/scripts/plot_splicing_results.R

echo "Generating Figure 4 (NMD Analysis)..."
run_r analysis/barcodeseq/rbp_barcode_screens/scripts/plot_nmd_results.R

echo "Generating Figure 5 (HHT Analysis)..."
run_r analysis/barcodeseq/rbp_barcode_screens/scripts/plot_eyfp_deopt_harr_results.R
run_r analysis/rnaseq/scripts/analyze_fold_changes.R deseq2_1.38.0
run_r analysis/riboseq/scripts/analyze_transcriptome_coverage.R

# Extended data figures
echo "Generating Extended Data figures..."
run_r analysis/barcodeseq/rbp_barcode_screens/scripts/compare_barcode_partitions.R
run_r analysis/flow_cytometry/figs1a_eyfp_reporter_sgeyfp_u2os_293t/scripts/plot_figs1a_flow.R
run_r analysis/flow_cytometry/figs1b_integration_efficiency_u2os_293t/scripts/plot_figs1b_flow.R
run_r analysis/polysome_profiling/figs2e_polysome_relic_hits/scripts/plot_figs2_polysomes.R
run_r analysis/qpcr/figs4_nmd_reporter_validation/scripts/plot_figs4_qpcr.R
run_r analysis/qpcr/figs5c_u937_gcn1_hht/scripts/plot_figs5c_qpcr.R
run_r analysis/qpcr/figs5d_zaki_gcn1_hht/scripts/plot_fig_s5d_qpcr.R

echo "
=== ANALYSIS COMPLETE ==="
echo "✓ All workflows executed successfully!"
echo "✓ All figures generated and saved to analysis/*/figures/"
echo "✓ Source data exported to source_data/"
echo ""
echo "See README.md for complete figure-to-code mapping"
echo "All done! 🎉"