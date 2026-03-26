#! /bin/bash

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
arg1 <- args[1]
arg1base <- basename(arg1)

#separate last column from averaging bedgraph files into a end and score column
oligomean <- read_tsv(arg1, col_names=c("chr", "start", "end"))
oligomean <- separate(oligomean, end, into =c("end", "score"), sep="\\|")
oligomean <- oligomean %>% mutate(end = as.double(end))
oligomean <- oligomean %>% mutate(score = as.double(score))

#remove 1bp from the end column so that there are no overlapping regions for bedgraph plotting in plotgardener
oligomean$end <- oligomean$end - 1

#write cleaned up mean bedgraph file
write_tsv(oligomean, paste0(arg1base, ".cleaned.bed"), col_names=FALSE)
