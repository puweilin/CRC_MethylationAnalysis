#The aim of this function is to plot the ROC curve for each of the variable as well as the combining all the variables with logistic regression. 
#Input1: @data
#Input2: @info: which is optional 
#Input3: @Selected.Genes: which of the variables should be plotted
#Input4: @CpGsite2Gene_method: which specifies the transformation method.



library(ggsci)
library(ggplot2)
library(ggthemes)
library(readxl)
library(reshape2)
library(readr)
library(dplyr)
library(tidyr)
library(ggbeeswarm)
library(gridExtra)
library(pROC)

q =  theme_classic() + 
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))+
  theme(legend.position=c(0.1,0.90))+
  theme(legend.title = element_blank())+
  theme(legend.text = element_text(colour="black", size = 16, face = "bold"))+
  theme(axis.text.x=element_text(size=15))+
  theme(axis.title.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  theme(axis.title.y=element_text(size=15))+
  theme(plot.title = element_text(hjust = 0.5, vjust=0))


ROCcurve = function(data, info=NULL, 
                    Selected.Genes, CpGsite2Gene_method, return.value=FALSE, saveplot=TRUE){
  
  new.data = CpGsite2Gene_method(data)
  data = new.data[[1]]
  Total.Genes = unique(data$Gene)
  data = subset(data, Gene %in% Selected.Genes)
  
  Sens = c()
  Spec = c()
  Gene = c()
  AUC = c()
  Figures = list()
  output_data = list()
  
  #get the roc results for multiple genes (n > 1)
  if(length(Selected.Genes) > 1){
    Selected.Genes = append(Selected.Genes, "Combined")
  }
  
  for(i in 1: length(Selected.Genes)){
    if(Selected.Genes[i] != "Combined"){
    temp.data = subset( data, Gene  == Selected.Genes[i] )} else{
      temp.data = data
    }
    Type = ifelse(substr(colnames(temp.data)[-c(1:3)],1,1) == "N", 0, 1)
    data_table = as.data.frame(t(temp.data[,-c(1:3)]))
    data_table = data.frame(Type, data_table)
    glm.fit = glm(Type ~. , data = data_table, family ="binomial")
    logistic.rocobj  = roc(na.omit(data_table)$Type, na.omit(as.vector(t(predict(glm.fit, data_table)))),
                           smooth = F)
    
    Sens = logistic.rocobj$sensitivities
    Spec = logistic.rocobj$specificities
    Gene = rep(Selected.Genes[i], length(logistic.rocobj$sensitivities))
    AUC =  format( logistic.rocobj$auc[[1]], digits = 2)
    print(sprintf("The AUC for the gene %s is %s", Selected.Genes[i], AUC))
    ROC.data = data.frame(Gene, Sens, Spec)
    output_data[[i]] = ROC.data
    
    if(saveplot){
      

    Best.idx = which.max(Sens + Spec)[1]
    Best.Sens = format(ROC.data$Sens[Best.idx], digits=2)
    Best.Spec = format(ROC.data$Spec[Best.idx], digits=2)
    Figures[[i]] = ggplot(ROC.data) + geom_path(aes(1-Spec, Sens), colour= pal_npg()(9)[1], lwd=3) +
      geom_abline(slope = 1, intercept = 0, lwd=3, colour="grey") + 
      q + scale_color_npg() + 
      ggtitle(Selected.Genes[i]) + theme(plot.title = element_text(hjust =0.5,size =20))+
      scale_y_continuous(limits = c(0,1.05), breaks = seq(0,1,0.2), expand = c(0.002,0.002)) +
      scale_x_continuous(limits = c(0,1.05), breaks = seq(0,1,0.2), expand = c(0.002,0.002)) +
      theme(legend.position = c(0.8, 0.2), legend.text = element_text(face="bold", size=18))+
      annotate("text", x = 0.7, y = 0.3, label = paste("AUC = ", AUC, sep=""), size=12)+
      annotate("text", x = 0.7, y = 0.2, label = paste("Sens = ", Best.Sens, sep=""), size=12)+
      annotate("text", x = 0.7, y = 0.1, label = paste("Spec = ", Best.Spec, sep=""), size=12)
    ggsave(Figures[[i]], filename = paste(Selected.Genes[i], "_ROC", ".pdf", sep=""), units = "cm", width = 18, height = 18,dpi = 300,device = "pdf")
    
    }
  }
  
  if(return.value){
    return(output_data)
  }else{
    ;
  }
  
 }