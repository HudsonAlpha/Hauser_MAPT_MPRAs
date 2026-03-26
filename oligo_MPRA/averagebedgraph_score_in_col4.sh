#! /bin/bash
#SBATCH -p normal
#SBATCH --mem=20GB
#SBATCH --job-name average_bedgraph
#SBATCH -o average_bedgraph%j.out

#This script takes a bedgraph with overlapping regions, and get the average score for each region
#This script assumes that you have your score in your 4th column. The tools require it to be in the 5th column
#if your file already has it in the 5th column, remove lines 19, 20, 30, 31, 34, and change 37 inputfile to $inputbed

#run by doing:
#sbatch averagebedgraph.sh bedgraph.bed

#output will be a sorted bedgraph file, a partitioned bed file, and a mean bedgraph file


inputbed=$1
basename=$(basename "$inputbed" .bed)
col5bed="${basename}.col5bed.bed"
col5bedintermediate="${basename}.col5.intermediate.bed"
sortedbed="${basename}.sorted.bed"
partitionedbed="${basename}.partition.bed"
meanbed="${basename}.sorted.mean.bed"

conda activate /cluster/home/rhauser/miniconda3/envs/myenv

####if your bedfile has the Score (-log10(pvalue)) in the 4th column.####

#bedmap expects it to be in the 5th column. This will move your score to column 5 and fill column 4 with NA
awk 'BEGIN {FS=OFS="\t"} {$5="NA"; print}' $inputbed > $col5bedintermediate
awk 'BEGIN {FS=OFS="\t"} {print $1, $2, $3, $5, $4}' $col5bedintermediate > $col5bed

#remove col5bedintermediate file
rm $col5bedintermediate

#file needs to be sorted first in order to work
sort -k1,1 -k2,2n $col5bed > $sortedbed
bedops -p $sortedbed > $partitionedbed
bedmap --echo --mean  $partitionedbed $sortedbed > $meanbed
