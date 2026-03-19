library(ggplot2)
library(tidyverse)
library(MPRAnalyze)


args <- commandArgs(trailingOnly = TRUE)
arg1 <- args[1]

dna_counts <- read.table("dna_counts_padded.csv", header = TRUE, sep=",", row.names=1)
rna_counts <- read.table("rna_counts_padded.csv", header = TRUE, sep=",", row.names=1)
dna_annot <- read.table("dna_annot.tsv", header = TRUE, sep="\t", row.names=1)
rna_annot <- read.table("rna_annot.tsv", header = TRUE, sep="\t", row.names=1)
mycondition <- arg1

#fill in NA with 0
dna_counts[is.na(dna_counts)] <- 0
rna_counts[is.na(rna_counts)] <- 0

#prep for filtering by DNA count of at least 50 per rep
#find columns for rep 1 and rep 2 DNA counts
dna1 <- grep("DNA_HEK293FT_1", colnames(dna_counts))
dna2 <- grep("DNA_HEK293FT_2", colnames(dna_counts))
dna3 <- grep("DNA_HEK293FT_3", colnames(dna_counts))

#calculate sums for rows in the identified columns
sumDNA1 <- rowSums(dna_counts[, dna1])
sumDNA2 <- rowSums(dna_counts[, dna2])
sumDNA3 <- rowSums(dna_counts[, dna3])

#filter rows for bins where both sums are greater than or equal to 50
filtered_rows <- (sumDNA1 >= 50) & (sumDNA2 >= 50) & (sumDNA3 >= 50)
filtered_dna_counts <- dna_counts[filtered_rows, ]
filtered_rna_counts <- rna_counts[filtered_rows, ]

#remove no_bc
filtered_dna_counts2 <- filtered_dna_counts[rownames(filtered_dna_counts) != "no_BC", ]
filtered_rna_counts2 <- filtered_rna_counts[rownames(filtered_rna_counts) != "no_BC", ]


myDnaCounts <- as.matrix(filtered_dna_counts2)
myRnaCounts <- as.matrix(filtered_rna_counts2)
myDnaAnnot <- dna_annot
myRnaAnnot <- rna_annot

myMPRAobject <- MpraObject(dnaCounts = myDnaCounts, rnaCounts = myRnaCounts, dnaAnnot = myDnaAnnot, rnaAnnot = myRnaAnnot)

myMPRAobject <- estimateDepthFactors(myMPRAobject, lib.factor = c("replicate", "condition"), which.lib="dna", depth.estimator ="uq")
myMPRAobject <- estimateDepthFactors(myMPRAobject, lib.factor = c("replicate", "condition"), which.lib="rna", depth.estimator ="uq")

myMPRAobject <- analyzeQuantification(obj = myMPRAobject, dnaDesign = ~ replicate)

alpha <- getAlpha(myMPRAobject, by.factor = "condition")

write.table(alpha, file="alpha_analysis2.tsv", quote=FALSE, sep='\t', col.names=NA)

png("alpha_boxplot_analysis2.png", width=800, height=800)
alphabox <- boxplot(alpha)
dev.off()


res.cond <- testEmpirical(obj = myMPRAobject, useControls = FALSE, statistic = alpha$V1, twoSided = FALSE)
summary(res.cond)
res.cond.labeled <- res.cond
res.cond.labeled$label <- row.names(alpha)

write.table(res.cond.labeled, file=paste0("testEmpirical_alpha_", mycondition, "_results_labeled_analysis2.tsv"), quote=FALSE, sep='\t', col.names=NA)

png("all_pval_hist_analysis2.png", width = 800, height = 800)
hist(res.cond$pval.mad, main="all")
dev.off()

