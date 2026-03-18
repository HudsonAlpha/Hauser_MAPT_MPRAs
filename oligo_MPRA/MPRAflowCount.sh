#! /bin/bash
#SBATCH -p normal
#SBATCH -o MPRAflow_count_%j.out

set +eu
source /cluster/home/rhauser/miniconda3/etc/profile.d/conda.sh
conda activate /cluster/home/rhauser/miniconda3/envs/MPRAflow

###### Edit These #####
MPRAflowDir="/cluster/home/rhauser/MPRAflow"
workingDir="/cluster/home/rhauser/files/rmhex02/MPRAflowCount_reps123/work_mpranalyze"
expFile="/cluster/home/rhauser/files/rmhex02/allfastqs/experiment_file_rmhex02_norep4.csv"
dataDir="/cluster/home/rhauser/files/rmhex02/allfastqs"
outputDir="/cluster/home/rhauser/files/rmhex02/MPRAflowCount_reps123/output_mpranalyze"
designFile="/cluster/home/rhauser/files/rmhex02/allfastqs/design_rmIllegalChars2.fa"
associationFile="/cluster/home/rhauser/files/rmhex02/allfastqs/RMHex02_bcassoc_allreads_filtered_coords_to_barcodes.pickle"
labelFile="/cluster/home/rhauser/files/rmhex02/allfastqs/rmhex02_label_rmIllegalChars_annotated.txt"
###########

cd $MPRAflowDir

# remove --mpranalyze flag if you want the basic mpraflow output instead

# I've been adding labels in R after the count script ran, but prior to mpranalyze,
# but to use the labels option you'll need to edit the labels file generated from
# the association script to give everything a group

nextflow run count_rmhedits.nf -w $workingDir --experiment-file $expFile --dir $dataDir --outdir $outputDir --design $designFile --association $associationFile --labels $labelFile --mpranalyze





