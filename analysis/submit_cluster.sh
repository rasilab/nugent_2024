snakemake \
    --jobs 999 \
    --cluster-config cluster.yaml \
    --cluster "sbatch -J {cluster.name} -n {cluster.cores}  -t {cluster.time} --mem {cluster.mem} -N {cluster.nodes}" \
    --use-singularity --latency-wait 10 \
    --singularity-args "--bind /hpc --bind /fh" \
    --cores=all \
    -p $@
