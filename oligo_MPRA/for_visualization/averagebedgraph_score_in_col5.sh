#! /bin/bash
#SBATCH -p normal
#SBATCH --mem=20GB
#SBATCH --job-name average_bedgraph
#SBATCH -o average_bedgraph%j.out

#This script takes a bedgraph with overlapping regions, and get the average score for each region
#This script assumes that you have your score in your 5th column. The tools require it to be in the 5th column

#run by doing:
#sbatch averagebedgraph.sh bedgraph.bed

#output will be a sorted bedgraph file, a partitioned bed file, and a mean bedgraph file


inputbed=$1
basename=$(basename "$inputbed" .bed)

sortedbed="${basename}.sorted.bed"
partitionedbed="${basename}.partition.bed"
meanbed="${basename}.sorted.mean.bed"

conda activate /cluster/home/rhauser/miniconda3/envs/myenv



#file needs to be sorted first in order to work
sort -k1,1 -k2,2n $inputbed > $sortedbed
bedops -p $sortedbed > $partitionedbed
bedmap --echo --mean  $partitionedbed $sortedbed > $meanbed
