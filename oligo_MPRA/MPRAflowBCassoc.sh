#!/bin/bash
#SBATCH -p normal
#SBATCH -o MPRAflow_assoc_%j.out

###change these variables###
MPRAflowDir="/cluster/home/rhauser/MPRAflow"
workingDir="/cluster/home/rhauser/files/rmhex02/MPRAflowBCassoc_piecebypiece_notstrandspecific/working"
outputDir="/cluster/home/rhauser/files/rmhex02/MPRAflowBCassoc_piecebypiece_notstrandspecific/output"
read1="/cluster/home/rhauser/files/rmhex02/barcode_association_Jane/mergeallfastqs_nodemux/allread1merge.fastq.gz"
read2="/cluster/home/rhauser/files/rmhex02/barcode_association_Jane/mergeallfastqs_nodemux/allread2merge.fastq.gz"
bcRead="/cluster/home/rhauser/files/rmhex02/barcode_association_Jane/mergeallfastqs_nodemux/allbarcodemerge.fastq.gz"
designFile="/cluster/home/rhauser/files/rmhex02/MPRAflowBCassoc_piecebypiece_notstrandspecific/design_rmIllegalChars_onelineseq.fa"
name="assoc_basic"

######

cd $MPRAflowDir

set +eu
source /cluster/home/rhauser/miniconda3/etc/profile.d/conda.sh
conda activate /cluster/home/rhauser/miniconda3/envs/MPRAflow

nextflow run association_rmhedits2.nf --w $workingDir --fastq-insert $read1 --fastq-insertPE $read2 --fastq-bc $bcRead --design $designFile --name $name --outdir $outputDir
