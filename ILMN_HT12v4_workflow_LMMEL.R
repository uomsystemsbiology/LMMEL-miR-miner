# setwd to contain sample/control probe data
setwd("/Volumes/Seagate Expansion Drive/R_Scripts_Repository/Scripts_for_miR-mRNA_paper/")
library(limma)
library(gplots)

## read in raw data ----
x <- read.ilmn(files="iltp3783_Sample_Probe_Profile.txt", ctrlfiles="Control Probe Profile.txt", other.columns="Detection")

## background correct, quantile normalize, log transform using negative and positive control probe information ----
y <- neqc(x, negctrl="NEGATIVE", regular="regular", offset=16) ##note default offset=16 adds numeric value of 16 to everything (ie: 4 after log2 transformation)

## read in E or M annotation and tidy data ----
samples <- read.delim("EorMordered.txt")
datanorm <- y[,-50] #remove SKMEL28 as not part of LM-MEL panel
groups <- t(samples[-50]) #three groups: E=epithelial M=mesenchymal C=control
colnames(groups) <- "phenotype"
colours <- groups
colours <- gsub("E", "red", colours)
colours <- gsub("M", "blue", colours)
colours <- gsub("C", "black", colours)
plotMDS(datanorm, labels=groups, col=colours)

## construct the linear models and perform DE analysis ----
design <- model.matrix(~0+groups)
fit <- lmFit(datanorm, design)
myContrasts <- makeContrasts(MvsE=groupsM-groupsE,
                              MvsC=groupsM-groupsC,
                              EvsC=groupsE-groupsC, levels=design)
fit2 <- contrasts.fit(fit, myContrasts)
fit3 <- eBayes(fit2)

## first contrast: M vs E ----
tT_MvsE <- topTable(fit3, coef=1, adjust.method="BH", p.value=0.05, sort.by="p", resort.by="logFC", number=Inf)
MvsE_expression <- datanorm$E[row.names(tT_MvsE),]

bluegreen <- colorRampPalette(c("deepskyblue1", "blue3", "green4", "green1"), bias=1)(100)
heatmap.2(MvsE_expression, dendrogram="column", scale="row", trace="none",
          margins=c(3,2), col=bluegreen,
          ColSideColors=(colours),
          key=TRUE, keysize=1.2, density.info="none", key.xlab="standardised intensity", key.ylab=NA,
          labRow=NA,
          srtCol=0, adjCol=c(0.5,NA), cexCol=1.2)
		  
MvsE_DEgeneList <- cbind(tT_MvsE, MvsE_expression)

## other contrasts as above, coef=2 or 3
