# The aim of this script is to filter the CpGsites and annotate these CpGsites with gene names. It is completed by using the Reference
# file named "SelectedRegion.xlsx", which tells the actual genomic coordiates of the targeted regions. And we want to exclude the CpGsites
# that is far from the targeted regions, which is uncorrectly mapped to the unwanted regions. 

#@Output: MethyTable.CG.RData
library(readxl)
library(readr)
RegionCompare = function(x, y,left, right, chr){
  if(x >left & x<right & y ==chr){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

MultiRegionCompare = function(x, y){
  POSstatus = apply(SelectedRegion, 1, function(a) RegionCompare(as.numeric(x), as.character(y), 
                                                                 as.numeric(a[3]),as.numeric(a[4]),as.character(a[2])))
  if(any(POSstatus==TRUE)){
    seq = which(POSstatus == TRUE)
    Genes = as.character(SelectedRegion[seq,1])
  }else{
    Genes = "NotIn"
  }
  return(Genes)
}


load("MethyTable.RData")
MethyTable.CG = subset(MethyTable, MethyTable$DiContext=="CG")
#Type = sapply(MethyTable.CG$Sample, function(x){ return(substr(as.character(x),1,1))})
#Type = ifelse(Type=="T", "Cancer","Control")

#select the regions that is close to the targeted genes
MethyTable.CG = subset(MethyTable.CG, MethyTable.CG$TotalReads >=15)
Num = sapply(MethyTable.CG$Sample, function(x) as.numeric(stri_extract( x, regex = "[0-9]+")) )
Type = ifelse(Num %% 2 == 0, "Control", "Cancer")
MethyTable.CG = data.frame(Type, MethyTable.CG)

SelectedRegion = read_excel("SelectedRegions.xlsx")
SelectedRegion$Gene = as.character(SelectedRegion$Gene)
SelectedRegion$Chr = as.character(SelectedRegion$Chr)

Gene = apply(MethyTable.CG, 1, function(b) MultiRegionCompare(as.numeric(b[5]),b[3]))
MethyTable.CG = data.frame(MethyTable.CG , Gene)
MethyTable.CG = subset(MethyTable.CG, MethyTable.CG$Gene != "NotIn")
save(MethyTable.CG, file="MethyTable.CG.RData")
