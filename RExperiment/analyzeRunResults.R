
require(ggiraph)
require(ggiraphExtra)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(corrplot)
df1=read.csv('/home/lu/AcrossTissue/csvs/5LS_phi_and_wormbaseId/merged_with_lower_bound/wholeGenomePhi_WBID.csv')
df2=read.csv('/home/lu/AcrossTissue/csvs/5LS_phi_and_wormbaseId/merged_with_lower_bound/5LS_empirical_exp.csv')

adult=read.csv('/home/lu/AcrossTissue/ROC_Runs/Run_L2L3_Combined/no_lower_bound_run_results/Organized/adult_Ce.fasta.wormbaseId_phi.csv')
dauer=read.csv('/home/lu/AcrossTissue/ROC_Runs/Run_L2L3_Combined/no_lower_bound_run_results/Organized/dauer_larva_Ce.fasta.wormbaseId_phi.csv')
emb=read.csv('/home/lu/AcrossTissue/ROC_Runs/Run_L2L3_Combined/no_lower_bound_run_results/Organized/elongating_embryo_Ce.fasta.wormbaseId_phi.csv')
L1=read.csv('/home/lu/AcrossTissue/ROC_Runs/Run_L2L3_Combined/no_lower_bound_run_results/Organized/L1_larva_Ce.fasta.wormbaseId_phi.csv')
L2L3=read.csv('/home/lu/AcrossTissue/ROC_Runs/Run_L2L3_Combined/no_lower_bound_run_results/Organized/L2L3_larva.fasta.wormbaseId_phi.csv')

df_selection=read.csv('/home/lu/AcrossTissue/ROC_Runs/Run_L2L3_Combined/no_lower_bound_run_results/Organized/selection.csv')
M <- cor(df_selection[,unlist(lapply(df_selection, is.numeric))  ])
corrplot(M, method = "number")




df= merge(df1,df2,by="WBID",all = TRUE)
df=merge(df,adult,by="WBID",all = TRUE)
df=merge(df,dauer,by="WBID",all = TRUE)
df=merge(df,emb,by="WBID",all = TRUE)
df=merge(df,L1,by="WBID",all = TRUE)
df=merge(df,L2L3,by="WBID",all = TRUE)
df
df=df[!duplicated(df$WBID), ]
df=df[!is.na(df$WBID), ]
df=df[!is.na(df$adult.Ce), ]
df=df[!is.na(df$Phi_WholeGenome), ]
summary(df)

#log
df_log=df
isnum <- sapply(df, is.numeric)
df_log[,isnum] <- lapply(df_log[,isnum], log1p)



M <- cor(df_selection[,unlist(lapply(df_selection, is.numeric))  ])
corrplot(M, method = "number")



summary(df_log)
summary(df)
colnames((df))
model <- lm(Phi_WholeGenome~elongating.embryo.Ce+L1.larva.Ce+dauer.larva.Ce+adult.Ce+L2L3_larva, data = df_log,na.action=na.omit)
model

model <- lm(Phi_WholeGenome~elongating.embryo.Ce+L1.larva.Ce+dauer.larva.Ce+adult.Ce+L2L3_larva, data = df)

length(df$Phi_WholeGenome)
model


length(dauer)

ggscatter(df, x="dauer.larva.Ce", y="PHI_Dauer") + yscale("log10", .format = TRUE)+xscale("log10", .format = TRUE)
ggscatter(df, x="L1.larva.Ce", y="PHI_L1") + yscale("log10", .format = TRUE)+xscale("log10", .format = TRUE)
ggscatter(df, x="adult.Ce", y="PHI_Adult") + yscale("log10", .format = TRUE)+xscale("log10", .format = TRUE)
ggscatter(df, x="L2L3_larva", y="PHI_L2L3") + yscale("log10", .format = TRUE)+xscale("log10", .format = TRUE)
ggscatter(df, x="elongating.embryo.Ce", y="PHI_Emb") + yscale("log10", .format = TRUE)+xscale("log10", .format = TRUE)

length(na.omit(df$PHI_Dauer))

#ggPredict(model,se=TRUE)


testData=df %>%select(elongating.embryo.Ce, L1.larva.Ce, dauer.larva.Ce,adult.Ce,L2L3_larva)
M_testData <- cor(testData)
corrplot(M_testData, method = "number")



pred.int <- predict(model,testData,se=T)
length(pred.int)
plot(df$Phi_WholeGenome,pred.int$fit)
fitValues=as.numeric(pred.int$fit)

cor(df$Phi_WholeGenome,df$adult.Ce,use="complete.obs",method="spearman")
cor(df$Phi_WholeGenome,df$dauer.larva.Ce,use="complete.obs",method="spearman")
cor(df$Phi_WholeGenome,df$L1.larva.Ce,use="complete.obs",method="spearman")
cor(df$Phi_WholeGenome,df$L2L3_larva,use="complete.obs",method="spearman")
cor(df$Phi_WholeGenome,df$elongating.embryo.Ce,use="complete.obs",method="spearman")
cor(df$Phi_WholeGenome,fitValues,use="complete.obs",method="spearman")


cor(df$PHI_Adult,df$adult.Ce,use="complete.obs",method="spearman")
cor(df$PHI_Dauer,df$dauer.larva.Ce,use="complete.obs",method="spearman")
cor(df$PHI_Emb,df$elongating.embryo.Ce,use="complete.obs",method="spearman")
cor(df$PHI_L1,df$L1.larva.Ce,use="complete.obs",method="spearman")
cor(df$PHI_L2L3,df$L2L3_larva,use="complete.obs",method="spearman")



#sort the df by "whole genome phi"
df_sorted <- df[order(-df$Phi_WholeGenome),]

emp_mean= rowMeans(subset(df_sorted, select = c(elongating.embryo.Ce,L1.larva.Ce,dauer.larva.Ce,adult.Ce,L2L3_larva)))

cor(emp_mean,df_sorted$Phi_WholeGenome,use="complete.obs",method="spearman")


