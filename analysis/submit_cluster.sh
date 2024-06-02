snakemake \
    --jobs 999 \
    --latency-wait 30 \
    --cluster-config cluster.yaml \
    --cluster "sbatch -n {cluster.n}  -t {cluster.time}" \
    --use-singularity \
    --use-envmodules \
    --singularity-args "--bind /fh --bind /hpc" \
    --cores=all \
    -p $@
