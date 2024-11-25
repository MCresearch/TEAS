#!/bin/bash
#SBATCH -J Si
#SBATCH -p gpu_4l
#SBATCH -N 1
#SBATCH -o nnof.out
#SBATCH -e nnof.err
#SBATCH --no-requeue
#SBATCH -A mhchen_g1
#SBATCH --qos=mhcheng4c
#SBATCH --gres=gpu:1
#SBATCH --overcommit
#SBATCH --mincpus=7

source /lustre3/mhchen_pkuhpc/mhchen_coe/4_dengzichao/env.sh
#/home/mhchen_pkuhpc/mhchen_coe/lustre2/1_sunliang/2_software/abacus-mlkedf-2.0/abacus-develop/source/module_hamilt_pw/hamilt_ofdft/ml_tools/build-noload/nnof > log
#/home/mhchen_pkuhpc/mhchen_coe/lustre2/1_sunliang/2_software/abacus-mlkedf-r_min/abacus-develop/source/module_hamilt_pw/hamilt_ofdft/ml_tools/build/nnof > log
export OMP_NUM_THREADS=1
python ../get_data.py
