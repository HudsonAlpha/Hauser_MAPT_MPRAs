#! /bin/bash
#SBATCH -p normal
#SBATCH -o MPRAflow_assoc_212_called_shorttime2_mapqcigar_%j.out


###change these variables###
workingDir="/cluster/home/rhauser/files/rmhex03/bcassoc/bcassoc_reseq/bcassoc212/workingcalledmapq"
outputDir="/cluster/home/rhauser/files/rmhex03/bcassoc/bcassoc_reseq/bcassoc212/outputcalledmapq"
designFile="/cluster/home/rhauser/files/rmhex03/design/varmpradesign_212bp.fasta"
name="assoc_basic_212_called_shorttime2_mapq"
cigar="212M"
read1="/cluster/home/rhauser/files/rmhex03/bcassoc/seq/fastqs_reseq_trimmed/RMHex03_bcassoc/RMHex03_bcassoc_i741_S1_R1_001.fastq.gz"
read2="/cluster/home/rhauser/files/rmhex03/bcassoc/seq/fastqs_reseq_trimmed/RMHex03_bcassoc/RMHex03_bcassoc_i741_S1_R3_001.fastq.gz"
bcRead="/cluster/home/rhauser/files/rmhex03/bcassoc/seq/fastqs_reseq_trimmed/RMHex03_bcassoc/RMHex03_bcassoc_i741_S1_R2_001.fastq.gz"
MPRAflowDir="/cluster/home/rhauser/MPRAflow"

######

cd $MPRAflowDir

set +eu
export HOST=$(hostname)
source /cluster/home/rhauser/miniconda3/etc/profile.d/conda.sh
conda activate /cluster/home/rhauser/miniconda3/envs/MPRAflow

nextflow run association_rmhedits2_shorttime2.nf --mapq 0 --cigar $cigar --w $workingDir --fastq-insert $read1 --fastq-insertPE $read2 --fastq-bc $bcRead --design $designFile --name $name --outdir $outputDir











