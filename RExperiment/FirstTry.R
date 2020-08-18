

library(SummarizedExperiment)
load("/home/lu/AcrossTissue/RExperiment/E-MTAB-2812-atlasExperimentSummary.Rdata")

rse <- experimentSummary$rnaseq



# write.table(fdf, "mydata.csv", sep=",")

#this gives the count information
countMatrix=assays(rse)$counts

# write.table(countMatrix, "/home/lu/AcrossTissue/RExperiment/count.csv", sep=",")

#This shows the col information
colData(rse)
#similar to the colData command
rse#dex
#a nice command to see the col names and row names
dimnames(rse)

rse#dex

srse=subset(rse, select = (sex == "hermaphrodite" & organism_part == "organism" & (genotype == "wild type genotype" |  genotype == "daf-2(e1370)III.")  ) )


# #create the metadata
# write.table(colData(srse), "/home/lu/AcrossTissue/RExperiment/metadata.csv", sep=",")
# #create the countMatrix
countMatrix=assays(srse)$counts
# write.table(countMatrix, "/home/lu/AcrossTissue/RExperiment/count.csv", sep=",")
library("DESeq2")

colData(srse)
library("airway")
data("airway")
se <- airway

countMatrix




srseDegObj=DESeqDataSet(rse,design = ~ 1)




dds <- DESeq(srseDegObj)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
rld <- rlog(dds)
plotPCA(rld)

# res <- results(dds)
# plotMA(res) #simple plot function
# # moderated log2 fold changes
# resultsNames(dds)
# resLFC <- lfcShrink(dds, coef=2, type="apeglm")
# # an alternate analysis: likelihood ratio test
# ddsLRT <- DESeq(dds, test="LRT", reduced= ~ 1)
# resLRT <- results(ddsLRT)

