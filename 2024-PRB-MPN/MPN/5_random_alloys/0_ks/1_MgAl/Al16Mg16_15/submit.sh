#!/bin/bash
#SBATCH -J Al16Mg16_15
#SBATCH -p cn_nl
#SBATCH -N 1 
#SBATCH -o relax.out
#SBATCH -e relax.err
#SBATCH --no-requeue
#SBATCH -A mhchen_g1 
#SBATCH --qos=mhchencnnl
#SBATCH -n 16

source /home/mhchen_pkuhpc/mhchen_coe/lustre2/1_sunliang/env.sh
export OMP_NUM_THREADS=1
mpirun -n 16 abacus > ../log/Al16Mg16_15 && echo "Al16Mg16_15 done"
