##########################################################################################
## This script is used to process the microRNA sequencing results from the LM-MEL panel ##
## following miRanalyzer processing. The mature miR read data is tabulated in files     ##
## called "mature_unique" and need to be placed in one folder (the wd) prior to running ##
## the following script. This script will output two summary datatables; one of the raw ##
## read counts and one of the RPKM data.                                                ##
##########################################################################################

## set up the function to extract the data from the data files

#extract the miRNA_ID and RPM data from an individual isoform quantification reads file
extractdata <- function(x, filepath1=filepath, datafiles1=datafiles,...) {
  reads <- read.delim(paste(filepath1, "/", datafiles1[x], sep=""), header=TRUE)
  reads <- reads[,c(1,3,5)]
  return(reads)
}

# setwd appropriately
setwd("#INSERT PATH TO DIRECTORY CONTAINING MATURE_UNIQUE DATA FILES#")

# make a list of all the files containing miRNA count data
filepath <- "./"
filelist <- dir(filepath)
datafiles <- filelist[grep("MEL", filelist)]

# process each file in "datafiles" to extract either read count (raw) or RPKM data
# append the results (the "summary" matrix) elementwise in a list
output_raw_list <- list(NA) #create a seed list
output_rpkm_list <- list(NA)

for (m in (1:length(datafiles))) {
  reads <- extractdata(m)
  reads_raw <- reads[,c(1,2)]
  output_raw_list[[m]] <- reads_raw
  reads_rpkm <- reads[,c(1,3)]
  output_rpkm_list[[m]] <- reads_rpkm
  rm(reads)
  rm(reads_raw)
  rm(reads_rpkm)
}
rm(m)

## now create a vector of the filenames (trimmed) to name each sample
sample_names <- gsub(".txt","", datafiles)

## then merge into summary dataframes for all detected miRNAs in all samples
raw_count_summary <- output_raw_list[[1]]
for (k in 2:length(datafiles)) {
  raw_count_summary <- merge(raw_count_summary, output_raw_list[[k]], by="name", all=TRUE)
} #ignore the warnings, they're just about non-unique column names whilst merging
rm(k)
colnames(raw_count_summary) <- c("miRNA_ID", sample_names)
write.table(raw_count_summary, file="#INSERT DESIRED FILENAME HERE#.txt", sep="\t")

rpkm_summary <- output_rpkm_list[[1]]
for (j in 2:length(datafiles)) {
  rpkm_summary <- merge(rpkm_summary, output_rpkm_list[[j]], by="name", all=TRUE)
} #ignore the warnings, they're just about non-unique column names whilst merging
rm(j)
colnames(rpkm_summary) <- c("miRNA_ID", sample_names)
write.table(rpkm_summary, file="#INSERT DESIRED FILENAME HERE#.txt", sep="\t")
