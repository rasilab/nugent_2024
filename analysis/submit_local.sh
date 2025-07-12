snakemake \
    --jobs 999 \
    --use-singularity \
    --singularity-args "--bind $(pwd)/../.." \
    --cores=all \
    -p $@
