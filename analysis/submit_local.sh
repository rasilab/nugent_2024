# Get the project root directory
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")"/../ && pwd)"
CONTAINER_DIR="${CONTAINER_DIR:-$PROJECT_ROOT/.env/singularity_cache}"

snakemake \
    --jobs 999 \
    --use-singularity \
    --singularity-args "--bind $(pwd)/../../.." \
    --config container_dir="$CONTAINER_DIR" \
    --cores=all \
    -p $@
