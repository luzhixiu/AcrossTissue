

library(SummarizedExperiment)
load("/home/lu/AcrossTissue/RExperiment/E-MTAB-2812-atlasExperimentSummary.Rdata")
rse <- experimentSummary$rnaseq
dim(rse)
assays(rse)
frse  <-  subset(rse, select = (sex == "hermaphrodite" & organism_part == "organism" & (genotype == "wild type genotype" |  genotype == "daf-2(e1370)III.")  ) )
unique(frse$sex)
fdf  <- assay(frse)
dim(fpd)
fdf  <- assay(frse)

write.table(fdf, "mydata.csv", sep=",")

#this gives the count information
countMatrix=assays(rse)$counts

write.table(countMatrix, "count.csv", sep=",")

