  GNU nano 5.6.1                                                                                                 rmhex01_count_call_shellscript_mpranalyze.sh                                                                                                            
#! /bin/bash
#SBATCH -p normal
#SBATCH -o MPRAflow_count_mpranalyze_%j.out

set +eu
source /cluster/home/rhauser/miniconda3/etc/profile.d/conda.sh
conda activate /cluster/home/rhauser/miniconda3/envs/MPRAflow

###### Edit These #####
MPRAflowDir="/cluster/home/rhauser/MPRAflow"
workingDir="/cluster/home/rhauser/files/rmhex01/MPRAflowCount_3reps/work_mpranalyze"
expFile="/cluster/home/rhauser/files/rmhex01/allrepsfastqs/experiment_file_rmhex01_20240909.csv"
dataDir="/cluster/home/rhauser/files/rmhex01/allrepsfastqs"
outputDir="/cluster/home/rhauser/files/rmhex01/MPRAflowCount_3reps/output_mpranalyze"
designFile="/cluster/home/rhauser/files/rmhex01/allrepsfastqs/100bindesignfile3.fa"
associationFile="/cluster/home/rhauser/files/rmhex01/allrepsfastqs/dictionary8.pickle"
###########

cd $MPRAflowDir

# remove --mpranalyze flag if you want the basic mpraflow output instead

# I've been adding labels in R after the count script ran, but prior to mpranalyze,
# but to use the labels option you'll need to edit the labels file generated from
# the association script to give everything a group

nextflow run count_rmhedits.nf -w $workingDir --experiment-file $expFile --dir $dataDir --outdir $outputDir --design $designFile --association $associationFile --mpranalyze













