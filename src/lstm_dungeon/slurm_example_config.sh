#!/bin/bash
#
#SBATCH -o /hpc/path/to/logs/out/slurm_%A_%a.out
#SBATCH -e /hpc/path/to/logs/err/slurm_%A_%a.err
#SBATCH --mem=48G
#SBATCH -p your-gpu-partition
#SBATCH --gres=gpu:1
#SBATCH --exclusive
#SBATCH -n 8
#SBATCH --array=1-30
#
echo ''
echo 'starting'
date
module load CUDA/11.4 #e.g.
srun /hpc/path/to/your/miniconda3/envs/nh2/bin/python /hpc/path/to/wd/lstms.py $SLURM_ARRAY_TASK_ID
echo 'done'
date

