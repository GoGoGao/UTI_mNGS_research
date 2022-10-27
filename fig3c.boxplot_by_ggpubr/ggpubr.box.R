library(ggplot2)
library(ggpubr)
library(patchwork)
library(dplyr)

file<-"ggpubr.box.input.table.xls"
dataset <- read.table(file,header = T,sep="\t")

dataset[dataset == 0] <-NA

for (i in 2:length(colnames(dataset))){
  df <- dataset[c(1,i)]
  df1 <-na.omit(df)
  my_comparisons <- list(c("Infected","Uninfected"))
  m = as.numeric(max(df1[,2]))
  v1 = 1.08*m

  p<- ggboxplot(df1, 
                x = "group_name", 
                y = colnames(dataset)[i],      
                color = "group_name", 
                palette = "npg" ,#两两比较的p值
                add = "jitter",#添加图形元素
				shape = "group_name"
  )+stat_compare_means(comparisons = my_comparisons,
                       label.y = c(v1)
  )+stat_compare_means(label.y = 1.5*m)+scale_colour_manual(values=c("#DC143C", "#20B2AA"))+
  theme(axis.text.x=element_text(vjust=1,size=14,color = "black"))+
  theme(axis.text.y=element_text(hjust=1,size=14,color = "black"))+
  theme(axis.title.y=element_text(hjust=0.5,size=14,color = "black"))+
  theme(axis.title.x=element_text(vjust=1,size=14,color = "black"))+
  xlab("")
  ggsave(p,file=paste(colnames(dataset)[i],".pdf",sep=""),width = 4,height =6,limitsize = F)
  ggsave(p,file=paste(colnames(dataset)[i],".png",sep=""),width = 4,height = 6,limitsize = F)
  write.table(df1,paste(colnames(dataset)[i],"filter.draw.input.xls",sep = ""),sep="	")
}

