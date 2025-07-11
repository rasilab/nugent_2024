snakemake \
    --jobs 999 \
    --use-singularity \
    --singularity-args "--bind /fh" \
    --cores=all \
    -p $@
