# This script is focused on combining the individual results into a complete methytable, which can be utilized for further analysis.
# It functions as merging the *.CGmap files and then do preliminary cleaning to remove the CpGsites that is not in the human chromosomes.
# And it saves the results as @MethyTable.RData
library(ggplot2)
library(ggthemes)
library(ggsci)
library(dplyr)
library(reshape)
library(readr)
library(readxl)


files.CGmap = list.files(pattern="*[.]CGmap")
Num = sapply(files.CGmap, function(x) as.numeric(strsplit(strsplit(x, '[.]')[[1]][1], 'A')[[1]][2]) )


Result.Merge=list()
for(i in 1:length(files.CGmap)){
  
  Result.Merge[[files.CGmap[i]]] =try(read.table(file=files.CGmap[i], sep="\t",header=F),silent=T)
  print(i)
}

Result = list()
for(i in 1:length(files.CGmap)){
  if(is.data.frame(Result.Merge[[i]]))
    Result[[files.CGmap[i]]] = Result.Merge[[i]]
}

for(i in 1:length(Result)){
  temp.data = Result[[i]]
  Sample = rep(strsplit(names(Result)[i],"[.]")[[1]][1], dim(temp.data)[1])
  temp.data = data.frame(Sample, temp.data)
  Result[[i]]= temp.data
}

MethyTable = bind_rows(Result)
colnames(MethyTable) = c("Sample", "CHR", "Strand", "POS", "Context", "DiContext","MethylRatio", "MethylReads", "TotalReads")
All.Chromos = c(seq(1,22,1),"X","Y","MT")
MethyTable = subset(MethyTable, MethyTable$CHR %in% All.Chromos)
save(MethyTable, file="MethyTable.RData")
