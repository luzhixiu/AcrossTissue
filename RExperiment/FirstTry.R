
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



countMatrix=assays(dds)$counts
write.table(countMatrix, "/home/lu/AcrossTissue/RExperiment/collasedReplicate.csv", sep=",",col.names=NA)
res <- results(dds)
dds$group <- dds$developmental_stage
dds$group <- factor(paste0(dds$developmental_stage,dds$AtlasAssayGroup))
dds <- DESeq(dds)


lfstages=colData(dds)$developmental_stage




ntd <- normTransform(dds)
library("vsn")
res <- results(ntd)


colData(ntd)
write.table(countMatrix, "/home/lu/AcrossTissue/RExperiment/log_normal_collasedReplicate.csv", sep=",",col.names=NA)

meanSdPlot(assay(ntd))

plotMA(res, ylim=c(-2,2))
library(ggfortify)


df <- read.csv(file = '/home/lu/AcrossTissue/RExperiment/collasedReplicate.csv',header = TRUE, row.names = 1, sep = ",")
??result

colData(dds)
results(dds,contrast = c("developmental_stage","X4.cell.embryo.Ce","newly.molted.young.adult.hermaphrodite.Ce"))


