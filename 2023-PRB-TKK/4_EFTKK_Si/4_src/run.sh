#!/bin/bash
#SBATCH -J EFTKK0
#SBATCH -p cn_nl
#SBATCH -N 1 
#SBATCH -o TKK.out
#SBATCH -e TKK.err
#SBATCH --no-requeue
#SBATCH -A mhchen_g1 
#SBATCH --qos=mhchencnnl
#SBATCH -n 24
source /appsnew/source/Python-3.7.3.sh
source /appsnew/source/openmpi-4.0.1-gcc.sh
python main.py
