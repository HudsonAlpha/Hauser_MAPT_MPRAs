#! /bin/bash

#SBATCH -p normal
#SBATCH --mem=35GB
#SBATCH -o vcfannot_%j.out

module load cluster/vep/112

FILE="$1"
OUT="annot_${FILE%.vcf}.tab"

vep --config newvep.ini \
    -i "$FILE" \
    -o "$OUT" \
    --force_overwrite \
    --tab \
    --fields "Location,Uploaded_variation,Allele,Existing_variation,CADD_PHRED"
