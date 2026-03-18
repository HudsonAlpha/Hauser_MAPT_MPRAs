  GNU nano 5.6.1                                                                                              bcassoc218_onlycalledreadsshorttime2_mapq0cigar.sh                                                                                                         
#! /bin/bash
#SBATCH -p normal
#SBATCH -o MPRAflow_assoc_218_called_shorttime2_%j.out


###change these variables###
MPRAflowDir="/cluster/home/rhauser/MPRAflow"
workingDir="/cluster/home/rhauser/files/rmhex03/bcassoc/bcassoc218/workingcalledmapq"
outputDir="/cluster/home/rhauser/files/rmhex03/bcassoc/bcassoc218/outputcalledmapq"
read1="/cluster/home/rhauser/files/rmhex03/bcassoc/seq/fastqs/RMHex03_bcassoc/RMHex03_bcassoc_i741_S1_R1_001.fastq.gz"
read2="/cluster/home/rhauser/files/rmhex03/bcassoc/seq/fastqs/RMHex03_bcassoc/RMHex03_bcassoc_i741_S1_R3_001.fastq.gz"
bcRead="/cluster/home/rhauser/files/rmhex03/bcassoc/seq/fastqs/RMHex03_bcassoc/RMHex03_bcassoc_i741_S1_R2_001.fastq.gz"
designFile="/cluster/home/rhauser/files/rmhex03/design/varmpradesign_218bp.fasta"
name="assoc_basic_218_called_shorttimemap1"
cigar="218M"

######

cd $MPRAflowDir

set +eu
export HOST=$(hostname)
source /cluster/home/rhauser/miniconda3/etc/profile.d/conda.sh
conda activate /cluster/home/rhauser/miniconda3/envs/MPRAflow

nextflow run association_rmhedits2_shorttime2.nf --mapq 0 --cigar 218M --w $workingDir --fastq-insert $read1 --fastq-insertPE $read2 --fastq-bc $bcRead --design $designFile --name $name --outdir $outputDir

