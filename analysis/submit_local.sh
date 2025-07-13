# Get the project root directory (handle symlinks properly)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR" && while [[ ! -f "run_everything.sh" ]]; do cd ..; done && pwd)"
CONTAINER_DIR="${CONTAINER_DIR:-$PROJECT_ROOT/.env/singularity_cache}"

echo "Using container directory: $CONTAINER_DIR"
snakemake \
    --jobs 999 \
    --use-apptainer \
    --apptainer-args "--bind $(pwd)/../../.." \
    --rerun-incomplete \
    --config container_dir="$CONTAINER_DIR" \
    --cores=all \
    -p $@
