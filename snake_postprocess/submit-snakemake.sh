#!/bin/bash
module unload python
module load Anaconda3/5.3.0
source activate /project2/nobrega/grace/conda/TWAS-env
module load python/2.7.12
export LD_LIBRARY_PATH=$HOME/bin/hdf5-1.8.20-linux-centos7-x86_64-gcc485-shared/lib/:$LD_LIBRARY_PATH

snakemake \
    --snakefile /project2/nobrega/grace/splicing/scripts/snake_postprocess/Snakefile \
    -kp \
    -j 30 \
    --rerun-incomplete \
    --cluster-config /project2/nobrega/grace/splicing/scripts/snake_postprocess/cluster.json \
    -c "sbatch \
        --mem={cluster.mem} \
        --nodes={cluster.n} \
        --tasks-per-node={cluster.tasks} \
        --partition=broadwl \
        --job-name={cluster.name} \
	   --output={cluster.logfile}" \
    $*
