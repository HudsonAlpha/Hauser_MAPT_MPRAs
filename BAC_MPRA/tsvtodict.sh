#! /bin/bash
#SBATCH -p normal
#SBATCH --mem=50GB

source /cluster/home/rhauser/miniconda3/etc/profile.d/conda.sh
conda activate /cluster/home/rhauser/miniconda3/envs/my_py_env

mydict=$1

python /cluster/home/rhauser/scripts/tsv_to_mpraflowdict_pickle.py $mydict
