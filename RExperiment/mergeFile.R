
require(ggiraph)
require(ggiraphExtra)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(corrplot)
df1=read.csv('/home/lu/AcrossTissue/csvs/5LS_phi_and_wormbaseId/merged_with_lower_bound/wholeGenomePhi_WBID.csv')
df2=read.csv('/home/lu/AcrossTissue/csvs/5LS_phi_and_wormbaseId/merged_with_lower_bound/5LS_empirical_exp.csv')

adult=read.csv('/home/lu/AcrossTissue/csvs/5LS_phi_and_wormbaseId/merged_with_lower_bound/adult_wormbaseId_phi.csv')
dauer=read.csv('/home/lu/AcrossTissue/csvs/5LS_phi_and_wormbaseId/merged_with_lower_bound/dauer_wormbaseId_phi.csv')
emb=read.csv('/home/lu/AcrossTissue/csvs/5LS_phi_and_wormbaseId/merged_with_lower_bound/emb_wormbaseId_phi.csv')
L1=read.csv('/home/lu/AcrossTissue/csvs/5LS_phi_and_wormbaseId/merged_with_lower_bound/L1_wormbaseId_phi.csv')
L2L3=read.csv('/home/lu/AcrossTissue/csvs/5LS_phi_and_wormbaseId/merged_with_lower_bound/L2L3_wormbaseId_phi.csv')

df_selection=read.csv('/home/lu/AcrossTissue/csvs/6LS_Selection.csv')
M <- cor(df_selection[,unlist(lapply(df_selection, is.numeric))  ])
corrplot(M, method = "number")


testData=df %>%select(elongating.embryo.Ce, L1.larva.Ce, dauer.larva.Ce,adult.Ce,L2L3_larva)
M_testData <- cor(testData)
corrplot(M_testData, method = "number")




df= merge(df1,df2,by="WBID",all = TRUE)
df=merge(df,adult,by="WBID",all = TRUE)
df=merge(df,dauer,by="WBID",all = TRUE)
df=merge(df,emb,by="WBID",all = TRUE)
df=merge(df,L1,by="WBID",all = TRUE)
df=merge(df,L2L3,by="WBID",all = TRUE)
df=df[!duplicated(df$WBID), ]
df=df[!is.na(df$WBID), ]
df=df[!is.na(df$adult.Ce), ]
df=df[!is.na(df$Phi_WholeGenome), ]
summary(df)

#log
df_log=df
isnum <- sapply(df, is.numeric)
df_log[,isnum] <- lapply(df_log[,isnum], log1p)

summary(df_log)
summary(df)
colnames((df))
model <- lm(Phi_WholeGenome~elongating.embryo.Ce+L1.larva.Ce+dauer.larva.Ce+adult.Ce+L2L3_larva, data = df_log,na.action=na.omit)
model




write.csv(df,"/home/lu/AcrossTissue/csvs/all_exp.csv")

write.csv(df_log,"/home/lu/AcrossTissue/csvs/all_exp_log.csv")

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
cor(df$Phi_WholeGenome,fitValues,use="complete.obs",method="spearman")

#sort the df by "whole genome phi"
df_sorted <- df[order(-df$Phi_WholeGenome),]

emp_mean= rowMeans(subset(df_sorted, select = c(elongating.embryo.Ce,L1.larva.Ce,dauer.larva.Ce,adult.Ce,L2L3_larva)))

cor(emp_mean,df_sorted$Phi_WholeGenome,use="complete.obs",method="spearman")

cor(emp_mean[1:100],df_sorted$Phi_WholeGenome[1:100],use="complete.obs",method="spearman")
cor(emp_mean[1:200],df_sorted$Phi_WholeGenome[1:200],use="complete.obs",method="spearman")
cor(emp_mean[1:300],df_sorted$Phi_WholeGenome[1:300],use="complete.obs",method="spearman")
cor(emp_mean[1:400],df_sorted$Phi_WholeGenome[1:400],use="complete.obs",method="spearman")
cor(emp_mean[1:500],df_sorted$Phi_WholeGenome[1:500],use="complete.obs",method="spearman")

















