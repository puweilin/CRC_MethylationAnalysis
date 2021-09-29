#The aim of this R code is to generate the Table 3 for the paper. The table 3 of the paper is designed to give a detailed 
# description of the prediction ability of the panel of several genes with some machine learning methods. For example, RF, SVM,
# NN, logistic regression, Naive Bayes, LDA, MDA, FDA, GBM, XGBoost, CatBoost etc. All the results given were 
# @Input: data
# @Input: info
# @Input: Selected.Genes
# @Input: table.name
# @Input: CpGsite2Gene_method
# @Input: N is for setting the times of computing with random training and validation datasets. 
# @Output: A Table

source("e:/R function/Methylation_Analysis_Pipeline//MLresult.R")
ML_classfiers_Table = function(data, info=NULL, Selected.Genes, table.name, CpGsite2Gene_method, N=1000){
   
  new.data = CpGsite2Gene_method(data,info)
  data = new.data[[1]]
  info = new.data[[2]]
  Total.Genes = unique(data$Gene)
  data = subset(data, Gene %in% Selected.Genes)
  Type = sapply(colnames(data)[-c(1:3)], function(x){ return(substr(x,1,1)) })
  Type = as.factor(ifelse(Type == "N", 0, 1))
  data_table = as.data.frame(t(data[,-c(1:3)])) 
  colnames(data_table) = Selected.Genes 
  methydata = data.frame(Type, data_table) 
  
  Sens.train=c()
  Spec.train=c()
  Accu.train=c()
  Sens.val  =c()
  Spec.val  =c()
  Accu.val  =c()
  require(caret)
  require(pROC)
  MLMatrix = matrix(nrow=N, ncol = 66)
  for(i in 1:N){
    methydata.case = subset(methydata, methydata$Type ==1)
    methydata.control = subset(methydata, Type ==0)
    train.idx.case = createDataPartition(methydata.case$Type, p =0.8, list=FALSE)
    train.idx.control= createDataPartition(methydata.control$Type, p =0.8, list=FALSE)
    trg.part = rbind( methydata.case[train.idx.case,], methydata.control[train.idx.control,])
    val.part = rbind( methydata.case[-train.idx.case,], methydata.control[-train.idx.control,])
    trg.part = na.omit(trg.part)
    val.part = na.omit(val.part)
    MLMatrix[i,] = MLresult(trg.part, val.part)
    if(i %% 100 ==0){
      print(i)
    }
  }
  MLMatrix = as.data.frame(MLMatrix)
  Name= c("Sens.log.tr", "Spec.log.tr", "Accu.log.tr", 
          "Sens.log.vl", "Spec.log.vl", "Accu.log.vl", 
          "Sens.rf.tr",  "Spec.rf.tr",  "Accu.rf.tr",  
          "Sens.rf.vl",  "Spec.rf.vl",  "Accu.rf.vl" ,
          "Sens.svm.tr", "Spec.svm.tr", "Accu.svm.tr", 
          "Sens.svm.vl", "Spec.svm.vl", "Accu.svm.vl",
          "Sens.NB.tr" , "Spec.NB.tr",  "Accu.NB.tr" , 
          "Sens.NB.vl" , "Spec.NB.vl" , "Accu.NB.vl" ,
          "Sens.NN.tr" , "Spec.NN.tr",  "Accu.NN.tr" ,
          "Sens.NN.vl" , "Spec.NN.vl",  "Accu.NN.vl" ,
          "Sens.LDA.tr", "Spec.LDA.tr", "Accu.LDA.tr",
          "Sens.LDA.vl", "Spec.LDA.vl", "Accu.LDA.vl",
          "Sens.MDA.tr", "Spec.MDA.tr", "Accu.MDA.tr",
          "Sens.MDA.vl", "Spec.MDA.vl", "Accu.MDA.vl",
          "Sens.FDA.tr", "Spec.FDA.tr", "Accu.FDA.tr",
          "Sens.FDA.vl", "Spec.FDA.vl", "Accu.FDA.vl",
          "Sens.GBM.tr", "Spec.GBM.tr", "Accu.GBM.tr",
          "Sens.GBM.vl", "Spec.GBM.vl", "Accu.GBM.vl",
          "Sens.xgb.tr", "Spec.xgb.tr", "Accu.xgb.tr",
          "Sens.xgb.vl", "Spec.xgb.vl", "Accu.xgb.vl",
          "Sens.catb.tr","Spec.catb.tr","Accu.catb.tr",
          "Sens.catb.vl","Spec.catb.vl","Accu.catb.vl")
  colnames(MLMatrix) = Name
  Avg.MLresult = as.vector(t(apply(MLMatrix, 2, mean)))
  Table = as.data.frame(matrix(Avg.MLresult,nrow=11, ncol=6, byrow=T))
  Table = format(Table, digits = 3)
  colnames(Table) =c("Train.Sens","Train.Spec","Train.Accu","Test.Sens","Test.Spec","Test.Accu")
  rownames(Table) =c("Logistic Regression","Random Forest", "SVM", "Naive Bayes","Neural Network",
                     "Linear Discriminant Analysis", "Mixture Discriminant Analysis", 
                     "Flexible Discriminant Analysis","Gradient Boosting Machine","XGBoost","CatBoost")
  #save(Table, file = paste(table.name, "RData", sep="."))
  write.csv(Table, paste(table.name, "csv", sep="."), row.names=T)
}
