# Author: Wanxin Wang
# Creates plots for 
# Fig.3A,B 
# Violin plots for Fig.3C, 5, S7

rm(list=ls())
library(Seurat)
library(RColorBrewer)
library(ggplot2)

immune_HQ_f<-readRDS("/data_share/immune_HQ_f.rds")
pal_ctp<-readRDS("/data_share/pal_ctp.rds")
nms.ctp_score<-readRDS("/data_share/nms.ctp_score.rds")

dirplot<-""

# Fig.3A UMAP
Idents(immune_HQ_f)<-"grps_high.res_68"
DimPlot(immune_HQ_f, cols = pal_ctp)

# Fig.3B DOTPLOT
DefaultAssay(immune_HQ_f)<-"RNA"
Idents(immune_HQ_f)<-"grps_high.res_68"
immune.markers_68<-FindAllMarkers(immune_HQ_f, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = T)

markers_high.res_68<-lapply(names(table(immune_HQ_f$grps_high.res_68)), function(x){
  cluster_markers<-immune.markers_68[immune.markers_68$cluster==x,]
  cluster_markers<-cluster_markers[order(cluster_markers$avg_logFC, decreasing = T),]
  #write.csv(cluster_markers, paste0(dir.data,"immune_clust",x,"_genes.csv"))
  cluster_markers
});names(markers_high.res_68)<-names(table(immune_HQ_f$grps_high.res_68))

# plot
ctp_genes_68<-unique(unlist(lapply(markers_high.res_68, function(x){x$gene[1:6]})))
gene_col<-rep("black",length(ctp_genes_68))
names(gene_col)<-ctp_genes_68
markers<-c("CD3D","CD3G","CD8A","CD8B","CD4","IGHM","MS4A1","CD160","KIT","CD14","NCAM1","FCGR3A")
markers<-markers[markers%in%ctp_genes_68]
gene_col[markers]<-"maroon"

pdf(paste0(dirplot,"dotplot_68.pdf"),width=9,height=3)
DefaultAssay(immune_HQ) <- "RNA"
plot<-DotPlot(immune_HQ_f, features = ctp_genes_68, group.by = "grps_high.res_68")
plot + theme(axis.text.y = element_text(angle = 0, size=8), axis.text.x = element_text(angle = 90, size=8, hjust=0.95,vjust=0.5, face="italic", colour = gene_col))+
  xlab("") + ylab("")&NoLegend()
dev.off()

# Fig.3C, 5CD, S7 VLNPLOTS
lapply(1:length(nms.ctp_score), function(x){
  nm<-nms.ctp_score[x]
  print(nm)
  f<-paste0(dirplot,"vlns/", nm, ".vln_68.pdf")
  pdf(f, width=3, height = 3.5)
  p<-VlnPlot(immune_HQ_f, group.by = "grps_high.res_68",cols=pal_ctp, features = nm, pt.size = 0)&
    NoLegend()&coord_flip()&geom_boxplot(width=0.1,fill="white",outlier.size = 0.4)&xlab("")&
    ylab("Signature score")&
    theme(axis.text.y=element_text(size=9,angle=0,hjust=1), 
          axis.text.x=element_text(size=9,vjust=1), 
          axis.title=element_text(size=10),
          plot.title = element_text(size=10))
  print(p)
  dev.off()
})
