# The aim of this R code is to automatically create to show the OR, 
# Confidence Interval, and P.value, Sensitivity, Specificity, AUC of each Gene with logistic regression.
# @Input: data
# @Input: Selected.Genes
# @Input: table.name
# @Input: CpGsite2Gene_method: which specifies the transformation method. 
# @Output: A Table
library(pROC)
MethylationStatus = function(data,  Selected.Genes, output_name, CpGsite2Gene_method){
  data = subset(data, Gene %in% Selected.Genes)
  new.data = CpGsite2Gene_method(data)
  data = new.data[[1]]
  data = as.data.frame(t(na.omit(t(data))),stringsAsFactors = F)
  data[,-c(1:3)] = apply(data[,-c(1:3)], 2, function(x) as.numeric(as.character(x)))
  
  Type = sapply(colnames(data)[-c(1:3)], function(x){ return(substr(x,1,1)) })
  Type = ifelse(Type == "N", 0, 1)
  data_table = as.data.frame(t(data[,-c(1:3)]))
  colnames(data_table) = data$Gene
  methydata = data.frame(Type, data_table)
  
  seq.case = which(methydata[,1] ==1)
  seq.control = which(methydata[,1] == 0)
  #Mean Case, Mean Control, Pvalue and Adjusted Pvalue
  if(dim(methydata)[2] >2){
  McaM = format( apply(methydata[,-1],2,function(x) {return( mean(x[seq.case], na.rm=TRUE))} ), digits = 2)
  McoM = format(apply(methydata[,-1],2,function(x) {return( mean(x[seq.control], na.rm=TRUE))} ), digits = 2)
  Pvalue=apply(methydata[,-1],2,function(x) {return( wilcox.test(x[seq.control], x[seq.case])$p.value)})
  FDR= format (p.adjust(Pvalue,method="fdr"), digits = 2, scientific = 2)
  Pvalue = format(Pvalue, digits =2, scientific  = 2)
  }else{
    McaM = format( mean(methydata[seq.case,-1], na.rm=TRUE), digits = 2)
    McoM = format( mean(methydata[seq.control,-1], na.rm=TRUE), digits = 2)
    Pvalue=wilcox.test(methydata[seq.control,-1], methydata[seq.case,-1])$p.value
    FDR = format (p.adjust(Pvalue,method="fdr"), digits = 2, scientific = 2)
    Pvalue = format(Pvalue, digits =2, scientific  = 2)
    
  }
  
  #Logistic regression analysis
  OR =c()
  CI.upper = c()
  CI.lower = c()
  Logistic.P = c()
  Sens=c()
  Spec=c()
  AUC =c()
  for(i in 1:(dim(methydata)[2] -1 )){
    temp = methydata[,c(1,i+1 )]
    temp[,1] = as.factor(temp[,1])
    glm.fit  = glm(temp[,1] ~ temp[,2], data = temp, family = "binomial") 
    OR[i] = format( log(exp(summary(glm.fit)$coefficients[2,1]),base = 10), digits =3 )
    Logistic.P[i] =format( summary(glm.fit)$coefficients[2,4], digits=2, scientific = 2)
    CI.upper[i] = format( log(exp(confint(glm.fit)[2,2]),base = 10), digits = 3)
    CI.lower[i] = format( log(exp(confint(glm.fit)[2,1]),base = 10), digits = 3)
    
    #Do the analysis of the sens, spec, and AUC
    predicted.value = predict(glm.fit)
    predicted.data  = data.frame(Type=na.omit(temp)[,1], predicted.value)
    logistic.rocobj  = roc(predicted.data$Type, predicted.data$predicted.value,smooth = FALSE)
    logistic.rocdata = data.frame(Sens = logistic.rocobj$sensitivities, Spec = logistic.rocobj$specificities)
    AUC[i] = format( logistic.rocobj$auc[[1]], digits = 2)
    #Find the best Sens and Spec
    logistic.rocdata[,3] = logistic.rocdata[,1] + logistic.rocdata[,2]
    seq.max = which(logistic.rocdata[,3] == max(logistic.rocdata[,3]))
    Sens[i] = format( logistic.rocdata[seq.max,1], digits = 2)
    Spec[i] = format( logistic.rocdata[seq.max,2], digits = 2)
    #print(i)
  }
  Logistic.P = format( p.adjust(Logistic.P, method = "fdr"), digits = 2, scientific = 2)
  Table = data.frame(Gene = data[,1], Chr = data[,2], Pos = data[,3], McaM, McoM, Pvalue, FDR, OR, CI.upper, CI.lower, Logistic.P, Sens,Spec, AUC)
  write.csv(Table, file = paste(output_name, "csv", sep="."), row.names=FALSE)
}
