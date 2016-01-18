## check that the microRNA and RNASeq data barcodes match

# setwd to parent directory containing untarred seq data
miRNAdata <- read.delim("./Output/TCGAset_mature_summary_barcodes.txt")
RNASeqdata <- read.delim("./Output/TCGAset_matched_RNASeq_summary_barcodes.txt")

miR_samples <- colnames(miRNAdata)[-1] #first column = miRNA_ID
RNASeq_samples <- colnames(RNASeqdata)[-1] #first column = gene_id

unique(miR_samples) #368
unique(RNASeq_samples) #368

match(miR_samples, RNASeq_samples)
# the patient-specific part is essentially encoded within the first three segments of the barcode
miR_samples_split <- strsplit(miR_samples, "\\.") #split the full barcode into segments separated by "."
# prepare a function to use with sapply on the above strsplit list
fixfunction <- function(x) {
  fixed <- paste(x[1], x[2], x[3], sep="-")
  return(fixed)
}
miR_samples_fixed <- sapply(miR_samples_split, FUN=fixfunction) #fix barcodes and return as a character vector

# do the same thing to the RNASeq barcodes
RNASeq_samples_split <- strsplit(RNASeq_samples, "\\.")
RNASeq_samples_fixed <- sapply(RNASeq_samples_split, FUN=fixfunction)
unique(miR_samples_fixed) #366
unique(RNASeq_samples_fixed) #366

# compare the two barcode vectors
matches <- match(miR_samples_fixed, RNASeq_samples_fixed)
matches[order(matches)]
# there are 2 patients who have two samples each in this dataset
counts <- matrix(data=NA, nrow=length(unique(miR_samples_fixed)), ncol=2) #dummy results matrix
for (i in 1:length(unique(miR_samples_fixed))) {
  counts[i,2] <- length(grep((unique(miR_samples_fixed))[i], miR_samples_fixed))
  counts[i,1] <- miR_samples_fixed[i]
}
counts[which(counts[,2]==2),1] #it's TCGA-ER-A2NF, and TCGA-GN-A4U7

## relabel the results dataframes and save them to text files
colnames(RNASeqdata) <- c("gene_id", RNASeq_samples_fixed)
write.table(RNASeqdata, file="./Output/RNASeq_summary_PatientBarcodes.txt", sep="\t")
colnames(miRNAdata) <- c("miRNA_id", miR_samples_fixed)
write.table(miRNAdata, file="./Output/miR_summary_PatientBarcodes.txt", sep="\t")