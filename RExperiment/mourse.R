
#This hunk of code makes the script run in multi-thread, doesn't affect the outcome of the code
library("BiocParallel")
register(MulticoreParam(4))

library(SummarizedExperiment)
load("/home/lu/AcrossTissue/RExperiment/E-MTAB-6798-mouse.Rdata")

rse <- experimentSummary$rnaseq


#this gives the count information
countMatrix=assays(rse)$counts

write.table(countMatrix, "/home/lu/AcrossTissue/RExperiment/mouse_count.csv", sep=",",col.names=NA)
colData(rse)
#This shows the col information
timeCollect=colData(rse)["age"]
#similar to the colData command
rse#dex
#a nice command to see the col names and row names
dimnames(rse)

rse#dex



#create the metadata
write.table(colData(rse), "/home/lu/AcrossTissue/RExperiment/mouse_metadata.SRR.csv", sep=",",col.names=NA)
#create the countMatrix
countMatrix=assays(rse)$counts
# write.table(countMatrix, "/home/lu/AcrossTissue/RExperiment/count.csv", sep=",",col.names=NA)
library("DESeq2")





srseDegObj=DESeqDataSet(rse,design = ~age)
dds <- DESeq(srseDegObj)
dds_bkup=dds
res <- results(dds)





#write the result to file
countMatrix=assays(dds)$counts
write.table(countMatrix, "/home/lu/AcrossTissue/RExperiment/collasedReplicate.csv", sep=",",col.names=NA)

































res
#filter out ones with low counts, genes with the sum of counts less than 10 are removed 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
rn=resultsNames(dds)
res <- results(dds)





#collapse the techinical replicates, grouped by the developmental stages
#dds <- collapseReplicates(dds, dds$developmental_stage,dds$technical_replicate_group)

res <- results(dds)
resultsNames(dds)

res <- results(dds, tidy=TRUE)


res <- res %>% mutate(sig=padj<0.05)

# How many of each?
res %>% 
  group_by(sig) %>% 
  summarize(n=n())
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="time")


plotDispEsts(dds)





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

summary(res)
lfstages=colData(dds)$developmental_stage




ntd <- normTransform(dds)

ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd),rank=F)
?meanSdPlot
