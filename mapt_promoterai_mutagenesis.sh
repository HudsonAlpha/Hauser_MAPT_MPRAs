#!/bin/bash

module load cluster/python/3.11.1
source promoterai/bin/activate

# chr17.fa is subset from the GRCh38 reference genome
promoterai --model_folder promoterAI_v1_hg38_mm10_finetune/ --var_file mapt_promoter_mutagenesis.tsv --fasta_file chr17.fa --input_length 20480

