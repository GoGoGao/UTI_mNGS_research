library(pheatmap)
x<-read.table("total.vfdb.even.rpm.sort.abseven.log10.xls",sep="\t",header=T,row.names=1)
group<-read.table("sort.group.list",sep="\t",header=F)
annotation_row = data.frame(Group=factor(group$V2))
rownames(annotation_row) = group$V1
ann_colors=list(Group=c("NC"="#808080","Infected"="#DC143C","Uninfected"="#20B2AA"))
pdf(file="heatmap.pdf",width=11,heigh=7)
pheatmap(t(x),fontsize=15,scale="none",cluster_cols=T,color = colorRampPalette(c("white","#9400D3"))(100),show_colnames=F,treeheight_row=80,treeheight_col=80,cluster_rows=F,annotation_row = annotation_row,annotation_colors=ann_colors,show_rownames=F,legend=T)   
dev.off()
