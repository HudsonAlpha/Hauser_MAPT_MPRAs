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
#filter counts for just mapt
maptcounts <- counts %>% filter(str_detect(name, regex("mapt", ignore_case = TRUE)))
maptcounts <- maptcounts %>% filter(!str_detect(name, regex("scram", ignore_case = TRUE)))

#load key
#key <- read_tsv("/cluster/home/rhauser/files/rmhex03/design/variant_reference_key.tsv")
#key$ALT <- key$id
#key <- key %>% dplyr::rename(ID = id)
#key <- key %>% dplyr::rename(REF = control)
#key <- key %>% select(ID, REF, ALT)

#load key
key <- read_tsv("/cluster/home/rhauser/files/rmhex03/design/fiveallelesmpradesign.tsv")

#filter key for just mapt
maptkey <- key %>% filter(str_detect(id, regex("mapt", ignore_case = TRUE)))
maptkey <- maptkey %>% filter(!str_detect(id, regex("scram", ignore_case = TRUE)))
newkey <- maptkey %>%
  group_by(element) %>%
  summarize(
    ID   = dplyr::first(element),
    REF  = id[alleleA],
        ALT1 = id[alleleB],
    ALT2 = id[alleleC],
    ALT3 = id[alleleD],
    ALT4 = id[alleleE],
    .groups = "drop"
  )


create_var_df_satmut <- function(df, map_df) {
  # Check required columns
  required_cols <- c("ID", "REF", "ALT1", "ALT2", "ALT3", "ALT4")
  if (!all(required_cols %in% colnames(map_df))) {
    stop("map_df must contain columns: 'ID', 'REF', 'ALT1', 'ALT2', 'ALT3', 'ALT4'")
  }

  if (!"name" %in% colnames(df)) {
    stop("df must contain column 'name'")
  }

  # Reshape map_df to long format for ALTs
  alt_cols <- c("ALT1", "ALT2", "ALT3", "ALT4")
  map_alt_long <- map_df %>%
    pivot_longer(cols = all_of(alt_cols), names_to = "alt_type", values_to = "ALT")

  # Check for matching refs or alts
  if (!any(df$name %in% map_df$REF) & !any(df$name %in% map_alt_long$ALT)) {
    stop("No matches found between the 'name' column in 'df' and the 'REF'/ALT columns in 'map_df'.")
  }

  # Merge on REF
  df_ref <- merge(df, map_df[, c("ID", "REF")], by.x = "name", by.y = "REF", all.x = FALSE)
  df_ref$allele <- "ref"

  # Merge on ALT (long format)
  df_alt <- merge(df, map_alt_long[, c("ID", "ALT", "alt_type")], by.x = "name", by.y = "ALT", all.x = FALSE)
  df_alt$allele <- df_alt$alt_type
  df_alt$alt_type <- NULL
  # Remove ALT from ref to avoid confusion
  df_ref$REF <- NULL

  # Combine results
  df_combined <- dplyr::bind_rows(df_ref, df_alt)

  # Select and return final output
  var_df <- df_combined %>%
    dplyr::select(variant_id = ID, allele, Barcode, dplyr::matches("count"))

  return(var_df)
}

vardf <- create_var_df_satmut(maptcounts, newkey)

dna_var <- create_dna_df(vardf)
rna_var <- create_rna_df(vardf)

mympra <- MPRASet(DNA = dna_var, RNA = rna_var, eid = row.names(dna_var), barcode = NULL)

nreps <- 4
bcs <- ncol(dna_var) / nreps
design <- data.frame(intcpt = 1, alt1 = grepl("ALT1", colnames(mympra)), alt2 = grepl("ALT2", colnames(mympra)), alt3 = grepl("ALT3", colnames(mympra)), alt4 = grepl("ALT4", colnames(mympra)))
block_vector <- rep(1:nreps, each=bcs)
mpralm_fit_var <- mpralm(object = mympra, design = design, aggregate = "none", normalize = TRUE, model_type = "corr_groups", plot = FALSE, block = block_vector)

top_var_1 <- topTable(mpralm_fit_var, coef = 2, number = Inf)
top_var_2 <- topTable(mpralm_fit_var, coef = 3, number = Inf)
top_var_3 <- topTable(mpralm_fit_var, coef = 4, number = Inf)
top_var_4 <- topTable(mpralm_fit_var, coef = 5, number = Inf)

top_var_1 <- top_var_1 %>%
  rownames_to_column(var="element")
top_var_2 <- top_var_2 %>%
  rownames_to_column(var="element")
top_var_3 <- top_var_3 %>%
  rownames_to_column(var="element")
top_var_4 <- top_var_4 %>%
  rownames_to_column(var="element")

top_var_1$allele <- "B"
top_var_2$allele <- "C"
top_var_3$allele <- "D"
top_var_4$allele <- "E"
toptab_mapt <- rbind(top_var_1, top_var_2, top_var_3, top_var_4)

#add element names back in


write_tsv(toptab_mapt, "bcalm_mapt_heks_10bcsthres.tsv")
