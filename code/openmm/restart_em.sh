#!/usr/bin/env bash
# Set walltime limit
#BSUB -W 72:00
#
# Set output file
#BSUB -o holo.nonoverlay.6.restarted.out
#
# Set error file
#BSUB -eo holo.nonoverlay.6.restarted.stderr 
#
# Specify node group
#BSUB -m "ls-gpu lt-gpu lp-gpu lg-gpu"
#BSUB -q gpuqueue
#
# nodes: number of nodes and GPU request
#BSUB -n 1 -R "rusage[mem=12]"
#BSUB -gpu "num=1:j_exclusive=yes:mode=shared"
#
# job name (default = name of script file)
#BSUB -J "holo nonoverlay 6 restarted"

python restart_em_state.py "/data/chodera/zhangi/vir_collaboration/data/em_output/holo.nonoverlay.6.50ns.minimized.pdb"
