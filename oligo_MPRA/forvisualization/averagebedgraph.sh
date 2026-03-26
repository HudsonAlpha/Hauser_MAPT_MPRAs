#! /bin/bash
#SBATCH -p normal
#SBATCH --mem=20GB
#SBATCH --job-name average_bedgraph
#SBATCH -o average_bedgraph%j.out

#This script takes a bedgraph with overlapping regions, and get the average score for each region

#run by doing:
#sbatch averagebedgraph.sh bedgraph.bed

#output will be a sorted bedgraph file, a partitioned bed file, and a mean bedgraph file

inputbed=$1
basename=$(basename "$inputbed" .bed)
sortedbed="${basename}.sorted.bed"
partitionedbed="${basename}.partition.bed"
meanbed="${basename}.sorted.mean.bed"

conda activate /cluster/home/rhauser/miniconda3/envs/myenv

#uncomment this if your bedfile has the Score (-log10(pvalue)) in the 4th column.
#bedmap expects it to be in the 5th column. This will move your score to column 5 and fill column 4 with NA
#this will overwrite your bed file, so save a spare in a different location just in case until this is done
#alternatively, you can run these lines before running this script to get a new bedgraph file with the score in column 5
#awk 'BEGIN {FS=OFS="\t"} {$5="NA"; print}' $inputbed > $inputbed
#awk 'BEGIN {FS=OFS="\t"} {print $1, $2, $3, $5, $4}' $inputbed > $inputbed

#file needs to be sorted first in order to work
sort -k1,1 -k2,2n $inputbed > $sortedbed
bedops -p $sortedbed > $partitionedbed
bedmap --echo --mean  $partitionedbed $sortedbed > $meanbed
