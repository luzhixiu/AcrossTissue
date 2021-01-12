



library(ggplot2)


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
colnames(df)

write.csv(df,"all_exp.csv")
