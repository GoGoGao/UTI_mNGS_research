args<-commandArgs(T)
args_infile = args[1]
species_infile = args[2]
pdfname = args[3]
procrustes_cluster1<-paste(pdfname,".procrustes.pdf",sep='')

#ARG
library(vegan)
arg_file <- read.table(args_infile,header = T,sep="\t",row.names = 1)
arg_file_trans <- as.data.frame(t(arg_file))
#Species
species_file <- read.table(species_infile,header = T,sep="\t",row.names = 1)
species_file_trans <- as.data.frame(t(species_file))

#dist
arg.dist <- vegdist(arg_file_trans)
species.dist <- vegdist(species_file_trans)

#降维
mds_arg <- monoMDS(arg.dist)
mds_species <- monoMDS(species.dist)

#pu
set.seed(1)
pro.s.r <- protest(mds_species,mds_arg,permutations=1999,symmetric=TRUE)
M2 <- sprintf("M^2 == %.4f",pro.s.r$ss)
pvalue <- sprintf("pvalue == %.4f",pro.s.r$signif)
#procru_test <- sprintf("italic(Procrustes analysis\:)")
labels <- data.frame(M2=M2,pvalue=pvalue,stringsAsFactors = FALSE)

library(ggplot2)

#获得X和y轴的坐标及旋转过的坐标
Pro_Y <- cbind(data.frame(pro.s.r$Yrot),data.frame(pro.s.r$X))
Pro_Y$group <- c(rep('Infected',30),rep('Uninfected',12))
Pro_Y$group <- as.factor(Pro_Y$group)

Pro_X <- data.frame(pro.s.r$rotation)

plot <- ggplot(data=Pro_Y) + geom_segment(aes(x=X1,y=X2,xend=MDS1,yend=MDS2),color = "gray",size=1.5)+
  geom_point(aes(X1,X2,shape = "ARG",color=group),se=FALSE,size=5)+
  geom_point(aes(MDS1,MDS2,shape = "Species",color=group),se=FALSE,size=5)+
  theme(panel.grid=element_blank(),panel.background = element_rect(color = 'black',fill = 'transparent'),legend.key = element_rect(fill='transparent'),
        axis.ticks.length = unit(0.4,"lines"),axis.ticks = element_line(color="black"),axis.line = element_line(color="black"),
        axis.title.x = element_text(colour = 'black',size = 16),axis.title.y = element_text(colour = 'black',size = 16),
        axis.text = element_text(colour = 'black',size = 14)) +
  labs(x='Dimesnsion 1',y='Dimesnsion 2',color = 'black')+
  geom_text(data=labels,mapping=aes(x = -0.15,y=0.17,label=M2),parse = TRUE,inherit.aes = FALSE,size = 6)+
  geom_text(data=labels,mapping=aes(x = -0.15,y=0.14,label=pvalue),parse = TRUE,inherit.aes = FALSE,size = 6)+
  theme(plot.title = element_text(size=14,colour="black",hjust=0.5,face="bold"))+  theme(legend.title = element_text(size=15),legend.text = element_text(size=10))

cairo_pdf(filename=procrustes_cluster1,height=10,width=12)
plot + scale_colour_manual(name="group",values = c("#00AFBB","#E7B800")) + scale_shape_manual(name="batch",values = c(16,17))
dev.off()






