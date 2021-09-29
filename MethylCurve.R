#The aim of this function is to plot the methylation curve plot based on the methylation percentage of the samples 
#Input1: @data: fortmat as shown in sample.xlsx, Tumor sample starts with T, Normal samples starts with N;
#Input2: @info: optional
#Input3: @subgroup: indicating if using part of the case samples for visualization;
#Input4: @subgroup_name: The name of the subgroup;
#Input5: @subgroup_samples: The samples included in the subgroup;
#Input6: @case.name: the case type name;
#Input7: @control.name: the control type name;
#Input8: @return.value: TRUE or FALSE, indicating that if the users want the figures as the return value.
#Input9: @legend.x: a vector, indicating the x axis of the legend for each plot
#Input10: @legend.y: a vector, indicating the y axis of the legend for each plot
#Input11: @error.bar: The type of error bar
library(dplyr)
library(tidyr)

q =  theme_classic() + 
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))+
  theme(legend.position=c(0.1,0.90))+
  theme(legend.title = element_blank())+
  theme(legend.text = element_text(colour="black", size = 16, face = "bold"))+
  theme(axis.text.x=element_text(size=15, angle = 20, vjust = 0.5))+
  theme(axis.title.x=element_text(size=20))+
  theme(axis.text.y=element_text(size=15))+
  theme(axis.title.y=element_text(size=20, vjust = 2))+
  theme(plot.title = element_text(hjust = 0.5, vjust=0))

MethylCurve = function(data, info=NULL, subgroup = FALSE, subgroup_name = "", subgroup_samples = list(),
                       control.name="", case.name="", return.value=FALSE, legend.x, legend.y, error.bar ="se"){
  
  genes = unique(data$Gene)
  
  Figures = list()
  for(i in 1:length(genes)){
    gene.idx = grep(genes[i], data$Gene)
    temp.data = data[gene.idx, ]
    temp.data = temp.data %>% gather( key = "Sample", value="Percentage", -Gene, -Chr, - Pos)
    Type = ifelse( sapply(temp.data$Sample, function(x) {substr(x,1,1)}) == "N", control.name, case.name)
    temp.data = data.frame(Type, temp.data)
    temp.data$Type = as.character(temp.data$Type)
    
    if(subgroup){
      for(j in 1:length(subgroup_name)){
        temp.idx = which(temp.data$Sample %in% subgroup_samples[[j]])
        temp.data$Type[temp.idx] = subgroup_name[j]
      }
      temp.data = subset(temp.data, Type %in% c(subgroup_name, control.name))
    }
    
    
    Mean = temp.data%>% group_by(Pos,Type) %>% dplyr::summarise(Mean = mean(Percentage, na.rm=TRUE))
    Sd =   temp.data%>% group_by(Pos,Type) %>% dplyr::summarise( Sd = sd(Percentage, na.rm=TRUE))
    N =  temp.data %>% dplyr::group_by(Pos, Type) %>% dplyr::summarise( N = dplyr::n())
    Upper.CI = temp.data %>% group_by(Pos, Type) %>% summarise(CI = t.test(Percentage, na.rm=T)$conf.int[2])
    Lower.CI = temp.data %>% group_by(Pos, Type) %>% summarise(CI = t.test(Percentage, na.rm=T)$conf.int[1])
    if(error.bar =="sd"){
      Sd = Sd$Sd
      CI1 = Mean$Mean-Sd
      CI2 = Mean$Mean+Sd
    }else if(error.bar == "se"){
      Se = Sd$Sd/sqrt(2* N$N)
      CI1 = Mean$Mean-Se
      CI2 = Mean$Mean+Se
    }else if(error.bar =="CI"){
      CI1 = Lower.CI$CI
      CI2 = Upper.CI$CI
    }else{
      print("Please select the right value for error bar!")
      break;
    }
    temp.plotdata = data.frame(Mean, CI1, CI2)
    colnames(temp.plotdata) = c("Position", "Type", "Mean", "CI1", "CI2")
    
    if(subgroup){
      temp.plotdata$Type = factor(temp.plotdata$Type, levels = c(subgroup_name, control.name))
    }else{
      temp.plotdata$Type = factor(temp.plotdata$Type, levels = c(case.name, control.name))
    }
    
    x.limits = c(min(temp.plotdata$Position)-20, max(temp.plotdata$Position)+20)
    Figures[[i]] = ggplot(temp.plotdata, aes(Position, Mean, colour=Type, group=Type)) + 
      geom_line(lwd = 1.5)+
      geom_point(size = 5) + 
      geom_errorbar(aes(ymin= CI1, ymax=CI2), size = 1, width=20,
                    position=position_dodge(0.05), colour = "black") + q +
      ylab("Mean Methylation Rate")+
      annotate("text", x= mean(x.limits), y = max(temp.plotdata$CI2)+0.03 , label = genes[i], size = 10) +
      scale_y_continuous(limits=c(0, max(temp.plotdata$CI2)+0.15),
                         breaks=seq(0,max(temp.plotdata$CI2)+0.1, 0.1))+
      scale_x_continuous(limits=x.limits, 
                         breaks=seq(min(temp.plotdata$Position), max(temp.plotdata$Position), round((x.limits[2] - x.limits[1])/5))) +
      scale_color_npg()+ theme(legend.position = c(legend.x[i], legend.y[i])) 
    ggsave(Figures[[i]], filename = paste(genes[i], "_curve",".pdf",sep=""), units = "cm", width = 18, 
           height = 12,dpi = 300,device = "pdf")
  }
  
  if(return.value){
    return(Figures)
  }else{
    ;
  }
  
  
}

