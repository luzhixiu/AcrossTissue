df=read.csv("/home/lu/AcrossTissue/csvs/6LS_Selection.csv")
df
str(df)

library(corrplot)
cor = cor(df[3:9])
corrplot.mixed(cor, lower.col = "black", number.cex = .7)


#
library(ggplot2)
ggplot(data=df, aes(x=dose, y=len, fill=Codon)) +
  geom_bar(stat="identity")
# Use position=position_dodge()

df[3]

hist(df[3])


ggplot(df, aes(x = Codon, y = L1)) +
  geom_point(aes(color = factor()))
