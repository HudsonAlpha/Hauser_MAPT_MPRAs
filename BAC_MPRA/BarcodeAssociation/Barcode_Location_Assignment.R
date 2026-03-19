#Barcode_Location_Assignment_slurm_2023April28.R


#module load cluster/R/4.1.0
#cd /cluster/home/bmoyers/for_becky/barcode_location_assignment_2023April17


args <- commandArgs(trailingOnly=T)
theSegment <- as.numeric(as.character(args[1]))




outDir <- "/cluster/home/bmoyers/for_becky/barcode_location_assignment_2023April17/brokenUp/"
dataDir <- "/cluster/home/bmoyers/for_becky/barcode_location_assignment_2023April17/"



################################################################################
#1) Read in the table I previously created.
#1a) It should be in the correct order already, but confirm...
#2) For each barcode, identify all of the regions to which it maps.
#3) If the regions all turn out to be identical:
#3a) Record the number of reads which support this association.
#3b) Save the association.
#4) If the regions are non-unique:
#4a) Save the associations to a problem file.
################################################################################


saveFile <- paste(dataDir, "bedout_barcodeAssociations.bed", sep="")
finalBed <- read.table(saveFile, header=F, sep="\t", stringsAsFactors=F)
colnames(finalBed) <- c("chr", "start", "end", "read", "barcode")

totalLines <- nrow(finalBed)

if(theSegment==1) {
  start <- 1
}
if(theSegment!=1) {
  startSegment <- theSegment-1
  startSegment <- startSegment/100
  start <- startSegment*totalLines
  start <- start+1
}

if(theSegment==100) {
  end <- totalLines
}

if(theSegment!=100) {
  endSegment <- theSegment/100
  end <- endSegment*totalLines
}

finalBed <- finalBed[start:end,]


currentBarcode <- finalBed[1,"barcode"]
currentSet <- c()

totalLines <- nrow(finalBed)


saveFile_good <- paste(outDir, "barcode_location_support_table_", theSegment, ".txt", sep="")
saveFile_problematic <- paste(outDir, "barcode_location_support_table_problematic_", theSegment, ".txt", sep="")
for (i in 1:nrow(finalBed)) {
  if(i%%10000==0) {print(paste(i, nrow(finalBed)))}
  if(finalBed[i,"barcode"]==currentBarcode) {
    currentSet <- rbind(currentSet, finalBed[i,])
  }
  if(finalBed[i,"barcode"]!=currentBarcode || i == totalLines) {
    if(nrow(rbind(currentSet))>1) {
      currentSet_mini <- currentSet[,c(1,2,3,5)]
      numObs <- nrow(currentSet)
      currentSet_mini <- unique(currentSet_mini)
      if(nrow(currentSet_mini)==1) {
        currentSet_mini$numObs <- numObs
        write.table(currentSet_mini, saveFile_problematic, row.names=F, col.names=F, sep="\t", quote=F, append=T)
      }
      if(nrow(currentSet_mini)>1) {
        write.table(currentSet_mini, saveFile_good, row.names=F, col.names=F, sep="\t", quote=F, append=T)
      }
    }
    if(nrow(currentSet)==1) {
      currentSet_mini <- currentSet[,c(1,2,3,5)]
      numObs <- nrow(currentSet)
      currentSet_mini$numObs <- numObs
      write.table(currentSet_mini, saveFile_good, row.names=F, col.names=F, sep="\t", quote=F, append=T)
    }
    currentBarcode <- finalBed[i,"barcode"]
    currentSet <- finalBed[i,]
  }
}

#Fin.
