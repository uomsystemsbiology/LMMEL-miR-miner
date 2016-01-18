## This script is used to process the microRNA sequencing results from matched TCGA SKCM (melanoma)
## data files using the "isoform quantification" files which detail which miR, which genomic location,
## and whether the read maps to the mature, star, stemloop, precursor or unknown part of the gene

# setwd to parent directory containing untarred smallRNAseq data
## make a list of all the files containing miRNA count data ----
filepath <- "./miRNASeq/BCGSC_IlluminaHiSeq_miRNASeq/Level_3"
filelist <- dir(filepath)
datafiles <- filelist[grep("isoform.quantification.txt", filelist)] #n=368

## function to extract miRNA_ID and RPKM data from an individual isoform quantification reads file ----
extractdata <- function(x, filepath1=filepath, datafiles1=datafiles,...) {
  reads <- read.delim(paste(filepath1, "/", datafiles1[x], sep=""), header=TRUE)
  reads <- reads[,-c(3,5)]
  return(reads)
}

## function to remove the entries that aren't for mature or star forms ----
remove_others <- function(y, ...) {
  to_remove_precursor <- grep("precursor", y$miRNA_region)
  to_remove_stemloop <- grep("stemloop", y$miRNA_region)
  to_remove_unannot <- grep("unannotated", y$miRNA_region)
  to_remove <- c(to_remove_precursor, to_remove_stemloop, to_remove_unannot)
  to_remove <- to_remove[order(to_remove)]
  
  reads_cut <- y[-to_remove,]
  region <- as.character(reads_cut$miRNA_region)
  region[grep("mature", region)] <- "mature"
  region[grep("star", region)] <- "star"
  reads_cut$miRNA_region <- region
  colnames(reads_cut) <- c("miRNA_ID", "coords", "RPKMmapped", "region")
  return(reads_cut)
}

## function to collapse all mature or star reads from single genomic locations ----
collapse <- function(z, ...) {
  mirs <- levels(z$miRNA_ID)
  summary <- matrix(data=NA, nrow=length(mirs), ncol=3)
  colnames(summary) <- c("miR_ID", "matureRPKM", "starRPKM")
  for (i in (1:length(mirs))) {
    summary[i,1] <- mirs[i]
    summary[i,2] <- sum(z[(z$miRNA_ID==mirs[i] & z$region=="mature"),3])
    summary[i,3] <- sum(z[(z$miRNA_ID==mirs[i] & z$region=="star"),3])
  }
  return(summary)
}

## process each file in "datafiles", appending results elementwise in a list ----
output_list <- list(NA) #create a seed list
for (m in (1:length(datafiles))) {
  reads <- extractdata(m)
  reads_cut <- remove_others(reads)
  summary <- collapse(reads_cut)
  output_list[[m]] <- summary
}

## create a vector of the TCGA patient barcodes to name each sample
sample_barcodes <- gsub(".isoform.quantification.txt","", datafiles)

## separate the mature and star data into separate lists then merge into summary dataframes
## for all detected miRNAs in all tumours
mature_only <- lapply(output_list, FUN=function(x) {x<-x[,-3]})
mature_summary <- mature_only[[1]]
for (k in 2:length(mature_only)) {
  mature_summary <- merge(mature_summary, mature_only[[k]], by="miR_ID", all=TRUE)
} #warnings about non-unique temporary column names whilst merging
colnames(mature_summary) <- c("miRNA_ID", sample_barcodes)
write.table(mature_summary, file="./Output/TCGAset_mature_summary_barcodes.txt", sep="\t")

star_only <- lapply(output_list, FUN=function(x) {x<-x[,-2]})
star_summary <- star_only[[1]]
for (k in 2:length(star_only)) {
  star_summary <- merge(star_summary, star_only[[k]], by="miR_ID", all=TRUE)
} #warnings are about non-unique column names whilst merging
colnames(star_summary) <- c("miRNA_ID", sample_barcodes)
write.table(star_summary, file="./Output/TCGAset_star_summary_barcodes.txt", sep="\t")