## This script is used to process the RNA sequencing results from all available TCGA melanoma
## data files using the "*.rsem.genes.normalized_results" files which detail the gene_id
## normalized_count from the RNAseq data (ie: collapsed isoform counts).

# setwd to parent directory containing untarred RNAseq data
# make a list of all the files containing miRNA count data ----
filepath <- "./RNASeqV2/UNC_IlluminaHiSeq_RNASeqV2/Level_3"
filelist <- dir(filepath)
datafiles <- filelist[grep("rsem.genes.normalized_results", filelist)] #n=368

## to extract the normalized RNASeq counts from an individual rsem.genes.normalized_results file ----
extractdata <- function(x, filepath1=filepath, datafiles1=datafiles,...) {
  reads <- read.delim(paste(filepath1, "/", datafiles1[x], sep=""), header=TRUE)
  return(reads)
}

# process each file in "datafiles", append results elementwise in a list called "output_list" ----
output_list <- list(NA) #create a seed list
for (m in (1:length(datafiles))) {
  reads <- extractdata(m)
  output_list[[m]] <- reads
}

## merge the output_list into a summary dataframe (column names = source filenames) ----
RNASeq_summary <- output_list[[1]] #seed the dataframe structure
for (k in 2:length(output_list)) {
  RNASeq_summary <- merge(RNASeq_summary, output_list[[k]], by="gene_id", all=TRUE)
} #warnings about non-unique temporary column names
colnames(RNASeq_summary) <- c("gene_id", datafiles)

## change the filenames to the patient barcodes using the separate file mapping ----
# load in the file map (comes with the downloaded data archive)
filemap <- read.delim("./FILE_SAMPLE_MAP.txt", header=TRUE)
mapping <- match(datafiles, filemap[,1])
barcodes <- as.vector(filemap[mapping,2])
RNASeq_summary_barcodes <- RNASeq_summary
colnames(RNASeq_summary_barcodes) <- c("gene_id", barcodes)

# save the final output as a file ----
write.table(RNASeq_summary_barcodes, 
            file="./Output/TCGAset_matched_RNASeq_summary_barcodes.txt", sep="\t")