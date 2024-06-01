snakemake \
    --cores=all \
    --use-singularity \
    --singularity-args "--bind /fh --bind /hpc" \
    -p $@
