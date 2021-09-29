#The aim of this code is to merge the several CpGsites in a single gene to get an overall estimation for each sample and each gene.
#Here, we used the mean value as the representatives.
library(dplyr)
CpGsite2Gene_byMean = function(x, y =NULL, allow_missing_rate_of_CpGsite = 0.20){
 
    #remove the CpGsites with high missing rate 
    missing.percent = apply( x[,-c(1:3)], 1, function(x) sum(is.na(x))/length(x) )
    seq.remove = ifelse( all(missing.percent <= allow_missing_rate_of_CpGsite), NA,
                         as.numeric(which(missing.percent > allow_missing_rate_of_CpGsite)))
    x.subset = x[setdiff(1:dim(x)[1], seq.remove),]
    if( is.null(y)){
      ;
    }else{
      y.subset = subset(y, Pos %in% x$Pos[setdiff(1:dim(x)[1], seq.remove)])
      
    }
    
  #Get the cleaned data after filtering the CpGsites with high missing rates.   
    
  N.CpGs = x[,1:3] %>% group_by(Gene) %>% summarize(N = n())
  N.selected.CpGs = x.subset[,1:3] %>% group_by(Gene) %>% summarize(N = n())
  Chr = x[,1:2] %>% group_by(Gene) %>% summarize(Chr = sample(Chr,size = 1))
  Start = x[,1:3] %>% group_by(Gene) %>% summarize( Pos = min(Pos))
  End = x[,1:3] %>% group_by(Gene) %>% summarize( Pos = max(Pos))
  new.x.anno= right_join(N.CpGs,N.selected.CpGs, by="Gene")
  N.removed.CpGs = new.x.anno[,2] - new.x.anno[,3]

   mean.data = list()
  for(i in 1:length(new.x.anno$Gene)){
    temp = subset(x.subset, Gene == new.x.anno$Gene[i])
    temp = t(temp[,-c(1:3)])
    mean.data[[i]] = as.data.frame(apply(temp, 1, function(x) mean(x, na.rm = T)))
  }
  mean.data.merge = as.data.frame(t(bind_cols(mean.data)))
  colnames(mean.data.merge) = colnames(x)[-c(1:3)]
  
  #Combine these datasets and get the result: new.x
  Pos = paste(Chr$Chr, paste(Start$Pos, End$Pos, sep="-"), sep=":")
  new.x = data.frame(Gene = as.character(new.x.anno$Gene), paste(N.selected.CpGs$N, N.CpGs$N, sep="/"), Pos, mean.data.merge)
  colnames(new.x)[1:3] = c("Gene","N.selected/N.CpGs", "Pos")
  if(is.null(y)){
    new.y = NULL
  }else{
    new.y = y.subset %>% group_by(Sample, Gene) %>% summarise(MethylReads = sum(MethylReads), TotalReads = sum(TotalReads))
    new.y$Percentage = new.y$MethylReads/ new.y$TotalReads
  }
  
  return(list(new.x, new.y))
}
