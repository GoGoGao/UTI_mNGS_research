suppressPackageStartupMessages(library("optparse"))
library(linkET)
library(corrplot)
library(vegan)
library(ggcor)
library(ggplot2)
library(dplyr)

option_list<-list(
    make_option(c("-u","--otu"), help="The otu matrix"),
    make_option(c("-f","--factor"),help="The factors matrix"),
    make_option(c("-n","--name"),help="The outfile name")
)
opt<-parse_args(OptionParser(usage="%prog [options] file\n",option_list=option_list))
if(is.null(opt$otu) && is.null(opt$factor)){
    cat("Use %prog -h for more help info\nThe author: yuanlijuan\n")
    quit("no")
}
otu<-opt$otu  #otu数据格式为otu\tgroup\tsample1\tsample2\t...
factor<-opt$factor  #factor数据格式为sample name\tfactor1\tfactor2\t...

#读取输入数据
spe<-read.table(otu,row.names = 1,header = T,sep = "\t",check.names = F)
factors<-read.table(factor,row.names = 1,header = T,sep = "\t",check.names = F)
dim(spe)
dim(factors)
#利用otu数据中的group列构造分类list,用于mantel test分析
group_list = split(1:dim(spe)[1], spe$group)
spe_pre<-select(spe,-group)
spe2<-as.data.frame(t(spe_pre))

#mantel test分析
mantel <- mantel_test(spe2, factors, 
                      spec.select = group_list) %>%  
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf), 
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),#定义Mantel的R值范围标签，便于出图 
         pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), 
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
#写出mantel test结果
write.table(mantel,"mantel_test.xls",sep="\t",row.names=F)

#绘制理化数据热图
#spearman相关性
p<-quickcor(factors, type = "upper",method = "spearman",show.diag=F)+
  geom_square()+#定义成方块状,geom_circle2()可定义为圆圈状
  #geom_couple(aes(colour = pd, size = rd), data = mantel, curvature = 0.1) +
  anno_link(aes(colour = pd, size = rd), data = mantel,curvature = 0.1)+#定义连线 
  scale_size_manual(values = c(0.5, 1.2, 2.5))+ 
  #scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +  #相关性颜色标注相反，选择下一行命令自行设置颜色
  scale_fill_gradientn(colours = c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#F7F7F7","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))+  #设置热图颜色
  scale_colour_manual(values = c("#1B9E77","grey","#D95F02"))+  #设置连线颜色
   theme(axis.text = element_text(size = 16),legend.title=element_text(size=15),legend.text=element_text(size=14))+
  guides(size = guide_legend(title = "Mantel's r",#定义图例 
                             order = 2), 
         colour = guide_legend(title = "Mantel's p",  
                               order = 3), 
         fill = guide_colorbar(title = "Spearman's r", order = 4))
ggsave(paste(opt$name,"-corr-spearman",".png", sep=""), p, width = 22, height = 18)
ggsave(paste(opt$name,"-corr-spearman",".pdf", sep=""), p, width = 22, height = 18)
#pearson相关性
p<-quickcor(factors, type = "upper",method = "pearson",show.diag=F)+
  geom_square()+#定义成方块状,geom_circle2()可定义为圆圈状
  #geom_couple(aes(colour = pd, size = rd), data = mantel, curvature = 0.1) +
  anno_link(aes(colour = pd, size = rd), data = mantel,curvature = 0.1)+#定义连线 
  scale_size_manual(values = c(0.5, 1.2, 2.5))+ 
  #scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +  #相关性颜色标注相反，选择下一行命令自行设置颜色
  scale_fill_gradientn(colours = c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#F7F7F7","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))+  #设置热图颜色
  scale_colour_manual(values = c("#D95F02","#1B9E77","grey"))+  #设置连线颜色
  guides(size = guide_legend(title = "Mantel's r",#定义图例 
                             order = 2), 
         colour = guide_legend(title = "Mantel's p",
                               order = 3),
         fill = guide_colorbar(title = "Pearson's r", order = 4))
ggsave(paste(opt$name,"-corr-pearson",".png", sep=""), p, width = 22, height = 18)
ggsave(paste(opt$name,"-corr-pearson",".pdf", sep=""), p, width = 22, height = 18)
