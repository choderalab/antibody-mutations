#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 72:00
#
# Set output file
#BSUB -o  apo.%I.out
#
# Set error file
#BSUB -eo apo.%I.stderr 
#
# Specify node group
#BSUB -m "ls-gpu lt-gpu lp-gpu lg-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=5]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "apo[9-16]"

python run_em.py "/home/zhangi/choderalab/vir_collaboration/data/em_input/apo/apo."${LSB_JOBINDEX}".in"
