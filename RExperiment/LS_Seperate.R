library("DESeq2")

df <- read.csv(file = '/home/lu/AcrossTissue/RExperiment/collasedReplicate.csv',header = TRUE, row.names = 1, sep = ",")



dds <- DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~ cond)