
#This hunk of code makes the script run in multi-thread, doesn't affect the outcome of the code
library("BiocParallel")
register(MulticoreParam(4))

library(SummarizedExperiment)
load("/home/lu/AcrossTissue/RExperiment/E-MTAB-2812-atlasExperimentSummary.Rdata")

rse <- experimentSummary$rnaseq



# write.table(fdf, "mydata.csv", sep=",")

#this gives the count information
countMatrix=assays(rse)$counts

# write.table(countMatrix, "/home/lu/AcrossTissue/RExperiment/count.csv", sep=",",col.names=NA)

#This shows the col information
colData(rse)
#similar to the colData command
rse#dex
#a nice command to see the col names and row names
dimnames(rse)

rse#dex

srse=subset(rse, select = (sex == "hermaphrodite" & organism_part == "organism" & (genotype == "wild type genotype" |  genotype == "daf-2(e1370)III.")  ) )


# #create the metadata
# write.table(colData(srse), "/home/lu/AcrossTissue/RExperiment/metadata.csv", sep=",",col.names=NA)
# #create the countMatrix
countMatrix=assays(srse)$counts
# write.table(countMatrix, "/home/lu/AcrossTissue/RExperiment/count.csv", sep=",",col.names=NA)
library("DESeq2")

colData(srse)
countMatrix
srseDegObj=DESeqDataSet(rse,design = ~ 1)
dds <- DESeq(srseDegObj)

#filter out ones with low counts, genes with the sum of counts less than 10 are removed 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
resultsNames(dds)
res <- results(dds)

#collapse the techinical replicates, grouped by the developmental stages
dds <- collapseReplicates(dds, dds$developmental_stage,dds$technical_replicate_group)

#filter out ones with low counts, genes with the sum of counts less than 10 are removed 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
resultsNames(dds)
res <- results(dds)



countMatrix=assays(ddsColl)$counts
write.table(countMatrix, "/home/lu/AcrossTissue/RExperiment/collasedReplicate.csv", sep=",",col.names=NA)
res <- results(ddsColl)
colData(ddsColl)
sum(res$padj < 0.1, na.rm=TRUE)
p=plotMA(res)


ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))

write.table(res, "/home/lu/AcrossTissue/RExperiment/DEGSEQ_filter_result.csv", sep=",",col.names=NA)


# dds <- dds[keep,]
# rld <- rlog(dds)
# plotPCA(rld)

# res <- results(dds)
# plotMA(res) #simple plot function
# # moderated log2 fold changes
# resultsNames(dds)
# resLFC <- lfcShrink(dds, coef=2, type="apeglm")
# # an alternate analysis: likelihood ratio test
# ddsLRT <- DESeq(dds, test="LRT", reduced= ~ 1)
# resLRT <- results(ddsLRT)

