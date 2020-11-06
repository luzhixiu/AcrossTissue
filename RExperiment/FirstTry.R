
#This hunk of code makes the script run in multi-thread, doesn't affect the outcome of the code
library("BiocParallel")
register(MulticoreParam(4))

library(SummarizedExperiment)
load("/home/lu/AcrossTissue/RExperiment/E-MTAB-2812-atlasExperimentSummary.Rdata")

rse <- experimentSummary$rnaseq
rse
df=read.csv("/home/lu/AcrossTissue/RExperiment/WBID_Coding.csv")
#Filter the gene names based on the cds gene names Mike provided
rse=subset(rse, rownames(rse) %in% df$WBID_Coding)
rse
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
# write.table(colData(srse), "/home/lu/AcrossTissue/RExperiment/metadata.SRR.csv", sep=",",col.names=NA)
# #create the countMatrix
countMatrix=assays(srse)$counts
# write.table(countMatrix, "/home/lu/AcrossTissue/RExperiment/count.csv", sep=",",col.names=NA)
library("DESeq2")

colData(srse)
countMatrix



srseDegObj=DESeqDataSet(srse,design = ~0+developmental_stage )
dds <- DESeq(srseDegObj)
dds_bkup=dds
resultsNames(dds)
dds
dds=dds_bkup

#Validate it was co






#collapse the techinical replicates, grouped by the developmental stages
dds <- collapseReplicates(dds, dds$developmental_stage,dds$technical_replicate_group)


#write the result to file
countMatrix=assays(dds)$counts
write.table(countMatrix, "/home/lu/AcrossTissue/RExperiment/worm_collasedReplicate.csv", sep=",",col.names=NA)









#collapse the techinical replicates, grouped by the developmental stages
#dds <- collapseReplicates(dds, dds$developmental_stage,dds$technical_replicate_group)






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




