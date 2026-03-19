library(tidyverse)
library(ggplot2)
library(devtools)
#devtools::install_github("kircherlab/BCalm")
library(BCalm)

#loadcounts
counts <- read_tsv("/cluster/home/rhauser/MPRAsnakeflow/results/experiments/rmhex03Count/reporter_experiment.barcode.HEKs.fromFile.default.min_oligo_threshold_10.tsv")

#remove counts for regions with no label
counts <- counts %>% filter(oligo_name != "no_BC")

counts <- counts %>% dplyr::rename(name = oligo_name)
counts <- counts %>% dplyr::rename(Barcode = barcode)
#counts <- counts %>% dplyr::rename(dna_count_1 = 'DNA 1 (condition HEKs, replicate 1)')
#counts <- counts %>% dplyr::rename(dna_count_2 = 'DNA 2 (condition HEKs, replicate 2)')
#counts <- counts %>% dplyr::rename(dna_count_3 = 'DNA 3 (condition HEKs, replicate 3)')
#counts <- counts %>% dplyr::rename(dna_count_4 = 'DNA 4 (condition HEKs, replicate 4)')
#counts <- counts %>% dplyr::rename(rna_count_1 = 'RNA 1 (condition HEKs, replicate 1)')
#counts <- counts %>% dplyr::rename(rna_count_2 = 'RNA 2 (condition HEKs, replicate 2)')
#counts <- counts %>% dplyr::rename(rna_count_3 = 'RNA 3 (condition HEKs, replicate 3)')
#counts <- counts %>% dplyr::rename(rna_count_4 = 'RNA 4 (condition HEKs, replicate 4)')

#counts <- counts %>% dplyr::select(barcode, name, dna_count_1, rna_count_1, dna_count_2, rna_count_2, dna_count_3, rna_count_3, dna_count_4, rna_count_4)


#get variant mapping
#filter counts for just snvs
snvcounts <- counts %>% filter(str_detect(name, regex("SNV", ignore_case = TRUE)))

#load key
key <- read_tsv("/cluster/home/rhauser/files/rmhex03/design/variant_reference_key.tsv")
key$ALT <- key$id
key <- key %>% dplyr::rename(ID = id)
key <- key %>% dplyr::rename(REF = control)
key <- key %>% select(ID, REF, ALT)


#filter key for just snvs
snvkey <- key %>% filter(str_detect(ID, regex("SNV", ignore_case = TRUE)))

vardf <- create_var_df(snvcounts, snvkey)

dna_var <- create_dna_df(vardf)
rna_var <- create_rna_df(vardf)

mympra <- MPRASet(DNA = dna_var, RNA = rna_var, eid = row.names(dna_var), barcode = NULL)

nreps <- 4
bcs <- ncol(dna_var) / nreps
design <- data.frame(intcpt = 1, alt = grepl("alt", colnames(mympra)))
block_vector <- rep(1:nreps, each=bcs)
mpralm_fit_var <- mpralm(object = mympra, design = design, aggregate = "none", normalize = TRUE, model_type = "corr_groups", plot = FALSE, block = block_vector)

top_var <- topTable(mpralm_fit_var, coef = 2, number = Inf)

top_var <- top_var %>% rownames_to_column(var="element")

write_tsv(top_var, "bcalm_snvs_labeled_heks_10bcthres.tsv")

