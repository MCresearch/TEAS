#!/bin/bash
#SBATCH -J bcc10
#SBATCH -p cn_nl
#SBATCH -N 1 
#SBATCH -o relax.out
#SBATCH -e relax.err
#SBATCH --no-requeue
#SBATCH -A mhchen_g1 
#SBATCH --qos=mhchencnnl
#SBATCH -n 8

source /lustre3/mhchen_pkuhpc/mhchen_coe/4_dengzichao/env.sh
export OMP_NUM_THREADS=1
mpirun -n 8 /home/mhchen_pkuhpc/mhchen_coe/lustre2/1_sunliang/2_software/abacus-3.7/abacus-develop/build/abacus > ../log/bcc10 && echo "bcc10 done"
