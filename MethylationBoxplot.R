#The aim of this function is to plot the Figure 1 based on the specific CpGsite2Gene method and the specific gene sets. 
#Input1: @data: fortmat as shown in sample.xlsx, Tumor sample starts with T, Normal samples starts with N;
#Input2: @info: optional
#Input3: @Selected.Genes
#Input4: @subgroup: indicating if using part of the case samples for visualization;
#Input5: @subgroup_name: The name of the subgroup;
#Input6: @subgroup_samples: The samples included in the subgroup;
#Input7: @case.name: the case type name 
#Input8: @control.name: the control type name
#Input9: @CpGsite2Gene_method: which specifies the transformation method. 
#Input10: @return.value: TRUE or FALSE, indicating that if the users want the figures as the return value.
library(ggsci)
library(ggplot2)
library(ggthemes)
library(readxl)
library(reshape2)
library(readr)
library(dplyr)
library(tidyverse)
library(ggbeeswarm)
library(gridExtra)

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


MethylationBoxplot = function(data, info=NULL, Selected.Genes,  subgroup = FALSE, subgroup_name = "", subgroup_samples = list(), 
                              case.name="",control.name="", CpGsite2Gene_method, return.value=FALSE,...){
  
  gene.idx = which(data$Gene %in% Selected.Genes)
  rawdata = data[gene.idx, ]
  
  transformed_data = CpGsite2Gene_method(rawdata,info,...)
  new_data = transformed_data[[1]]
  Type = sapply(colnames(new_data)[-c(1:3)], function(x){ return(substr(x,1,1)) })
  Type = ifelse(Type == "N", control.name, case.name)
  data_table = as.data.frame(t(new_data[,-c(1:3)]))
  data_table = data.frame(Type, Sample = as.vector(t(colnames(rawdata)[-c(1:3)])), data_table)
  colnames(data_table) = c("Type","Sample", as.character(new_data$Gene))
  data_plot = data_table %>% gather(-Type, -Sample, key="Gene", value="Methylation")  
  data_plot$Type = as.character(data_plot$Type)
  
  if(subgroup){
    for(j in 1:length(subgroup_name)){
      temp.idx = which(data_plot$Sample %in% subgroup_samples[[j]])
      data_plot$Type[temp.idx] = subgroup_name[j]
    }
    data_plot = subset(data_plot, Type %in% c(subgroup_name, control.name))
  }
  Figures = list()
  for(i in 1: length(Selected.Genes)){
    idx = grep(Selected.Genes[i], data_plot$Gene)
    temp = data_plot[idx, ]
    
    if(subgroup){
      temp$Type = factor(temp$Type, levels = c(subgroup_name, control.name))
    }else{
      temp$Type = factor(temp$Type, levels = c(case.name, control.name))
    }
    temp$Methylation = as.numeric(temp$Methylation)
    Figures[[i]] = ggplot(data = temp) + geom_boxplot(aes(Type, Methylation),lwd=0.5,fill="lightgrey", outlier.shape = NA) + 
                   geom_quasirandom(aes(Type, Methylation, colour= Type),size=3) + q + scale_color_npg() +  
                   theme(axis.text.x=element_text(size=15)) + ggtitle( Selected.Genes[i]) + 
                   theme(axis.title.x=element_blank())+ theme(legend.position="none") + 
                   scale_y_continuous(limits=c(0,1), breaks = seq(0,1,0.2))
    ggsave(Figures[[i]], filename = paste(Selected.Genes[i],"_Boxplot",".pdf",sep=""), units = "cm", width = 16, height = 12,dpi = 300,device = "pdf")
    }
  if(return.value){
  return(Figures)
    }else{
    ;
  }
  
}
