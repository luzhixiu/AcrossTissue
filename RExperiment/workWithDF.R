
library(tidyverse)



df <- read.csv(file = '/home/lu/AcrossTissue/RExperiment/collasedReplicate.csv',header = TRUE, row.names = 1, sep = ",",stringsAsFactors=FALSE)

plot(df$X3.fold.embryo.Ce, log="y", type='h', lwd=10, lend=2)

df$dauer.larva.Ce
df
ls=df$X3.fold.embryo.Ce

plts=list()
names(df)
for (name in names(df)){
  ls=df[[name]]
  ls=log(ls)
  plt=hist(ls)
  plt
  append(plts,plt)
}



