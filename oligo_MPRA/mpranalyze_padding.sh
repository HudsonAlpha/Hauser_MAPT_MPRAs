#! /bin/bash

#SBATCH -p normal
#SBATCH --mem=10G
#SBATCH -o padding_%j.out


source activate mpranalyze

dnacounts="dna_counts.tsv"
output="dna_counts_padded.csv"
my_directory="/cluster/home/rhauser/scripts"

python $my_directory/mpranalyze_padding.py $dnacounts $output

#column count check
awk -F',' '{print NF}' $dnacounts | sort -n | uniq -c

rnacounts="rna_counts.tsv"
output2="rna_counts_padded.csv"

python $my_directory/mpranalyze_padding.py $rnacounts $output2

#column count check
awk -F',' '{print NF}' $rnacounts | sort -n | uniq -c
