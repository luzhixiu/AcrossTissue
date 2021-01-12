






df1=read.csv("/home/lu/Microbiome_ML/ColoradoData/LeavesAndRoot/formated/Root/ENV2014Roots_2020_06_07.csv")
df2=read.csv("/home/lu/Microbiome_ML/ColoradoData/LeavesAndRoot/formated/Root/ITS_roots_prop_at_least50.csv")


df= merge(df1,df2,by="Sample")
write.csv(df,"/home/lu/Microbiome_ML/ColoradoData/LeavesAndRoot/formatted/merge.csv")
