#! /bin/bash

#SBATCH --array 1-88%10
#SBATCH --mem=50G
#SBATCH -n 9

#first activate conda env: conda activate /cluster/home/aanderson/miniconda3/envs/idemux
#submit with path to fastqs file as arguments

samplesheet=$1

line=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${samplesheet})
current_sample=$(echo $line | awk '{print $1}')

echo $current_sample
./RNA_pipeline.sh $current_sample 


