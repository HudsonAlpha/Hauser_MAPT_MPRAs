#!/bin/bash

#SBATCH --array 1-88%10
#SBATCH --mem=20G
#SBATCH -n 9
#SBATCH -o gencounts.out

#submit with path to fastqs file as arguments

module load cluster/htslib/1.9
module load cluster/star/2.7.5c
module load cluster/samtools/1.16.1
module load cluster/bbmap/38.42
module load cluster/cutadapt

R1=$1
R2=$(echo $R1 | sed 's/_1.fq/_2.fq/')

sample=$(echo $R1 | cut -d "/" -f 11)

echo "R1: ${R1}"
echo "R2: ${R2}"
echo "Sample: ${sample}"

# symlink
mkdir -p fastq/${sample}
ln -s ${R1} fastq/${sample}/R1_fastq.gz
ln -s ${R2} fastq/${sample}/R2_fastq.gz


#unzip
mkdir -p fastq_unzip/${sample}
gunzip -c fastq/${sample}/R1_fastq.gz >  fastq_unzip/${sample}/R1_fastq
gunzip -c fastq/${sample}/R2_fastq.gz >  fastq_unzip/${sample}/R2_fastq


# trim
mkdir -p trimming/${sample}
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -m 30 -o trimming/${sample}/R1_fastq -p trimming/${sample}/R2_fastq fastq_unzip/${sample}/R1_fastq fastq_unzip/${sample}/R2_fastq

# align
mkdir -p alignment/${sample}/
STAR --runThreadN 9  --genomeDir /cluster/home/aanderson/myers/BrainTF/RNA/STAR_GENCODE_v42/ --readFilesIn trimming/${sample}/R1_fastq trimming/${sample}/R2_fastq --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 -->

# index
samtools index  alignment/${sample}/"Aligned.sortedByCoord.out.bam"


# dedup
mkdir -p deduplicated/${sample}/
java -jar /cluster/software/picard-2.21.6/picard.jar MarkDuplicates \
      I=alignment/${sample}/Aligned.sortedByCoord.out.bam \
      O=deduplicated/${sample}/removed_duplicates.bam \
      REMOVE_SEQUENCING_DUPLICATES=true \
      M=deduplicated/${sample}/marked_dup_metrics.txt
samtools sort deduplicated/${sample}/removed_duplicates.bam > deduplicated/${sample}/removed_duplicates.srt.bam
samtools index deduplicated/${sample}/removed_duplicates.srt.bam


# counting unique alignments
echo "counting"
mkdir -p htseq_counting/${sample}/


#htseq-count -m intersection-nonempty -s yes -f bam -r pos -t gene_id deduplicated/${sample}/removed_duplicates.srt.bam htseq_counting/${sample}/genes.gtf > htseq_counting/${sample}/htseq.unique.ID.tsv

htseq-count -m intersection-nonempty -s no -f bam -r pos -t gene deduplicated/${sample}/removed_duplicates.srt.bam /cluster/home/aanderson/myers/BrainTF/RNA/STAR/gencode.v42.primary_assembly.annotation.gtf  > htseq_counting/${sample}/htseq.unique.nonstranded.tsv



# clean
rm fastq_unzip/${sample}/R1_fastq
rm fastq_unzip/${sample}/R2_fastq

