
require(ggiraph)
require(ggiraphExtra)
require(plyr)


library(ggplot2)
library(ggpubr)

df1=read.csv('/home/lu/AcrossTissue/csvs/5LS_phi_and_wormbaseId/merged/wholeGenomePhi_WBID.csv')
df2=read.csv('/home/lu/AcrossTissue/csvs/5LS_phi_and_wormbaseId/merged/5LS_empirical_exp.csv')

adult=read.csv('/home/lu/AcrossTissue/csvs/5LS_phi_and_wormbaseId/merged/adult_wormbaseId_phi.csv')
dauer=read.csv('/home/lu/AcrossTissue/csvs/5LS_phi_and_wormbaseId/merged/dauer_wormbaseId_phi.csv')
emb=read.csv('/home/lu/AcrossTissue/csvs/5LS_phi_and_wormbaseId/merged/emb_wormbaseId_phi.csv')
L1=read.csv('/home/lu/AcrossTissue/csvs/5LS_phi_and_wormbaseId/merged/L1_wormbaseId_phi.csv')
L2L3=read.csv('/home/lu/AcrossTissue/csvs/5LS_phi_and_wormbaseId/merged/L2L3_wormbaseId_phi.csv')


df= merge(df1,df2,by="WBID",all = TRUE)
df=merge(df,adult,by="WBID",all = TRUE)
df=merge(df,dauer,by="WBID",all = TRUE)
df=merge(df,emb,by="WBID",all = TRUE)
df=merge(df,L1,by="WBID",all = TRUE)
df=merge(df,L2L3,by="WBID",all = TRUE)
df

#log
df_log=df
isnum <- sapply(df, is.numeric)
df_log[,isnum] <- lapply(df_log[,isnum], log1p)
df_log
summary(df_log)
summary(df)
colnames((df))
model <- lm(Phi_WholeGenome~elongating.embryo.Ce+L1.larva.Ce+dauer.larva.Ce+adult.Ce+L2L3_larva, data = df_log,na.action=na.omit)
model

write.csv(df_log,"/home/lu/AcrossTissue/csvs/all_exp_log.csv")

model <- lm(Phi_WholeGenome~elongating.embryo.Ce+L1.larva.Ce+dauer.larva.Ce+adult.Ce+L2L3_larva, data = df)
length(pred.int)
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
library(dplyr)

testData=df %>%
  select(elongating.embryo.Ce, L1.larva.Ce, dauer.larva.Ce,adult.Ce,L2L3_larva)

pred.int <- predict(model,testData,se=T)
length(pred.int)
plot(df$Phi_WholeGenome,pred.int)




