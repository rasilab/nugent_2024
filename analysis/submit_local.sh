snakemake \
    --jobs 999 \
    --use-singularity \
    --singularity-args "--bind /hpc --bind /fh" \
    --cores=all \
    -p $@
