library(tidyverse)
library(tibble)

args <- commandArgs(trailingOnly = TRUE)
arg1 <- args[1]


#load data from sh scripts as tibble

dat <- as_tibble(read.table(arg1, col.names = c("ID", "UMI1", "Index1", "UMI2", "Index2"), header = FALSE))

#dat <- as_tibble(read.table("headhekdna.all_info.txt", col.names = c("ID", "UMI1", "Index1", "UMI2", "Index2"), header = FALSE))
dat$spacer <- "+"

#merge ID row columns into one to make downstream easier
dat <- dat %>% select(ID, UMI1, Index1, spacer, UMI2) %>%
  unite("IDmerge", ID:UMI1, sep =" ", remove = TRUE) %>%
  unite("IDmerge2", IDmerge:Index1, sep = " ", remove = TRUE)


dat$fauxqual <- "AAAAAAAAAAAAAAAA"


write.table(dat[,c("IDmerge2", "UMI2", "spacer", "fauxqual")], 
            file = "UMI.fastq", 
            sep = "\n", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)

