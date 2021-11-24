#* Author: Wanxin Wang
#* This script calculates the scores and plot
#* scores in Fig.2B
#* including: permutation test results, ratioNext, onTarget, and median normalized signature scores

rm(list=ls())
library(Seurat)
library(RColorBrewer)
library(plyr)
library(ggplot2)
library(dplyr)

# ----------------------
# perform utest and
# calculate fc
# for two groups of cells
# ----------------------
u.fc<-function(ct,names1,names2,mth.padj="BH",fc.dm, alternative){
  uout<-utest_wrapper(ct,names1,names2, alternative=alternative)
  print(dim(uout))
  uout$"p.adj"<-p.adjust(uout$pval,method=mth.padj)
  print(dim(uout))
  
  fc<-apply(ct,1,function(x){
    fc<-(mean(x[names1])+fc.dm)/(mean(x[names2])+fc.dm)
    log2.fc<-log2(fc)
  })
  list("uout"=uout,"fc"=fc)
}

# ----------------------
# WRAPPER for u.fc()
# for calculating uniquely
# expressing genes in a
# one to all manner
# ----------------------
u.fc_one2all<-function(ct,names,mth.padj="BH",fc.dm, alternative){
  out<-lapply(c(1:length(names)), function(x){
    #print(x)
    print(length(names[[x]]))
    print(length(unlist(names[-x])))
    u.fc(ct, names[[x]],unlist(names[-x]),mth.padj=mth.padj,fc.dm=fc.dm, alternative=alternative)})
  out
}


dir.plot<-""
dir.o<-"/datashare/"

## import xcell signatures and sc data
xcell_sc<-readRDS("GSE111976_ct_endo_10x.rds") # can be downloaded from GSE111976
xcell_sig_all_list<-readRDS("/data_share/xcell_sig_all_list.rds")

# For each xCell signature and each single cell type
# calculate signature scores
sig_ctp<-names(xcell_sig_all_list)
ct<-xcell_sc[["RNA"]]@data

ctps<-names(table(xcell_sc$ctp_to_plot))
names_ctp<-lapply(ctps, function(x){
  names(xcell_sc$ctp_to_plot)[xcell_sc$ctp_to_plot==x]
});names(names_ctp)<-ctps

# calculate raw sig scores
xcellsig_vs_ctp<-lapply(1:length(sig_ctp),function(x){
  
  sig.type<-sig_ctp[x]
  print(paste0("Signature: ", sig.type))
  
  sig<-xcell_sig_all_list[[sig.type]]
  print(paste0("Length of signature: ", length(sig)))
  
  sig<-sig[sig%in%rownames(ct)]
  print(paste0("Number of genes in ct: ", length(sig)))
  
  ct.sig<-ct[sig,]
  
  # one2all utest, p, log2(FC)
  # alternative is "greater")
  u.fc_out<-u.fc_one2all(ct.sig,names=names_ctp,mth.padj="BH",fc.dm=1E-02,alternative="greater")
  names(u.fc_out)<-names(names_ctp)
  
  u.fc_out
});names(xcellsig_vs_ctp)<-sig_ctp

# For each xCell signature and each single cell type
# Calculate % genes in signature with p.adj < 0.05, log2FC > 1
list_pct_sigfig<-lapply(1:length(xcellsig_vs_ctp),function(i){
  o<-sapply(1:length(xcellsig_vs_ctp[[i]]),function(j){
    out<-xcellsig_vs_ctp[[i]][[j]]
    genes_padj<-rownames(out$uout)[out$uout$p.adj<0.05]
    genes_fc<-names(out$fc)[out$fc>1]
    genes_sigfig<-intersect(genes_padj,genes_fc)
    length_sigfig<-length(genes_sigfig)
    pct_sigfig<-length_sigfig/nrow(out$uout)
  },USE.NAMES = T)
  o
});names(list_pct_sigfig)<-names(xcellsig_vs_ctp)

mat_pct_sigfig<-as.data.frame(matrix(unlist(list_pct_sigfig),nrow=length(xcell_sig_all_list),byrow=T))
rownames(mat_pct_sigfig)<-names(list_pct_sigfig)
colnames(mat_pct_sigfig)<-names(xcellsig_vs_ctp[[1]])

#* ADJUST scores for epi subtypes to avoid over-penalizing against each other
#* that is, eliminate the other two epi subtypes when calculating epi subtype of interest
# for all signatures
sig_ctp<-names(xcell_sig_all_list)
ct<-xcell_sc[["RNA"]]@data
names_epi_subtypes<-names_ctp[epi_ctps]
xcellsig_vs_epi.by.subtypes<-lapply(1:length(sig_ctp),function(x){
  # define ct
  sig.type<-sig_ctp[x];print(paste0("Signature: ",sig.type))
  sig<-xcell_sig_all_list[[sig.type]]
  print(paste0("Length: ",length(sig)))
  sig<-sig[sig%in%rownames(ct)]
  print(paste0("Length in ct: ",length(sig)))
  ct.sig<-ct[sig,]
  
  #one2all utest, p, log2(FC)
  #IMPORTANT: alternative is "greater"
  u.fc_out<-lapply(1:3,function(i){
    u.fc(ct.sig,names1=names_epi_subtypes[[i]],names2=names_non.epi,mth.padj="BH",fc.dm=1E-02,alternative="greater")
  })
  names(u.fc_out)<-names(names_epi_subtypes)
  u.fc_out
});names(xcellsig_vs_epi.by.subtypes)<-sig_ctp

list_pct_sigfig_epi.by.subtypes<-lapply(1:length(xcellsig_vs_epi.by.subtypes),function(i){
  o<-sapply(1:length(xcellsig_vs_epi.by.subtypes[[i]]),function(j){
    out<-xcellsig_vs_epi.by.subtypes[[i]][[j]]
    genes_padj<-rownames(out$uout)[out$uout$p.adj<0.05]
    genes_fc<-names(out$fc)[out$fc>1]
    genes_sigfig<-intersect(genes_padj,genes_fc)
    length_sigfig<-length(genes_sigfig)
    pct_sigfig<-length_sigfig/nrow(out$uout)
  },USE.NAMES = T)
  o
});names(list_pct_sigfig_epi.by.subtypes)<-names(xcellsig_vs_epi.by.subtypes)

mat_pct_sigfig_epi.by.subtypes<-as.data.frame(matrix(unlist(list_pct_sigfig_epi.by.subtypes),ncol=3,byrow=T))
colnames(mat_pct_sigfig_epi.by.subtypes)<-names(xcellsig_vs_epi.by.subtypes[[1]])
rownames(mat_pct_sigfig_epi.by.subtypes)<-names(list_pct_sigfig_epi.by.subtypes)

# combined adjusted score for epi with that of other single cell cell types
mat_pct_sigfig_manuscript<-cbind(mat_pct_sigfig[,-grep("epithelium",colnames(mat_pct_sigfig))],
                                 mat_pct_sigfig_epi.by.subtypes[rownames(mat_pct_sigfig),])

# Normalize by row median
mat_pct_sigfig_manuscript_norm.med<-t(apply(mat_pct_sigfig_manuscript, MARGIN = 1, FUN = function(x){
  if(median(x)!=0){
    return(x/median(x))
  }else{
    return(x/(median(x)+0.01))
  }
}))
#saveRDS(mat_pct_sigfig_manuscript_norm.med, paste0(dir.o, "signature.score_norm.med.rds"))

# map between signature and cell type
map_ctp_target<-matrix(rep(0,nrow(mat_pct_sigfig)*ncol(mat_pct_sigfig)), nrow=nrow(mat_pct_sigfig))
colnames(map_ctp_target)<-colnames(mat_pct_sigfig)
rownames(map_ctp_target)<-rownames(mat_pct_sigfig)

map_ctp_target[grep("Epitheli",rownames(map_ctp_target)),grep("epitheli",colnames(map_ctp_target))]<-1
map_ctp_target[grep("Fibroblasts",rownames(map_ctp_target)),grep("fibroblast",colnames(map_ctp_target))]<-1
map_ctp_target[grep("Endothelial",rownames(map_ctp_target)),grep("Endothelium",colnames(map_ctp_target))]<-1
map_ctp_target[grep("CD4+|CD8+|Th|Tregs|Tgd|NK|NKT",rownames(map_ctp_target)),grep("Lymphocyte",colnames(map_ctp_target))]<-1
map_ctp_target[grep("Macrophages|DC|Monocytes",rownames(map_ctp_target)),grep("Macrophage",colnames(map_ctp_target))]<-1
map_ctp_target[grep("Smooth muscle|MSC|Pericyte|muscle",rownames(map_ctp_target)),grep("Smooth muscle",colnames(map_ctp_target))]<-1

#saveRDS(map_ctp_target, paste0(dir.o, "map_ctp_target.rds"))

#ratioNext + Specificity
mat_pct_sigfig_manuscript_spec<-lapply(1:nrow(mat_pct_sigfig_manuscript_norm.med), function(x){
  ctp<-rownames(mat_pct_sigfig_manuscript_norm.med)[x]
  print(ctp)
  map<-map_ctp_target[ctp,]
  ordered<-sort(mat_pct_sigfig_manuscript_norm.med[x,], decreasing = T)
  a<-ordered[1]
  #print(a)
  if(sum(map)==0){
    print("No target")
    up<-a
    dn<-ordered[2]
    tar=0
  }else{
    print("Target exists")
    target<-names(map)[map==1]
    non_target<-names(map)[map==0]
    b<-sort(ordered[non_target], decreasing = T)[1]
    if(b==a){
      # if the top hit occurs in a non-target cell type
      print("off-target!")
      tar=1
      up<-a
      dn<-ordered[2]

    }else{
      print("on-target!")
      tar=2
      up<-a
      dn<-b
    }
  }
  spec=up/dn
  return(list("specificity"=spec, "tar"=tar))
});names(mat_pct_sigfig_manuscript_spec)<-rownames(mat_pct_sigfig_manuscript_norm.med)

# ratioNext
score_specificity<-sapply(mat_pct_sigfig_manuscript_spec, function(x){
  x$specificity
});names(score_specificity)<-names(mat_pct_sigfig_manuscript_spec)
# saveRDS(score_specificity, paste0(dir.o, "specificity.score_ratio.rds"))

# onTarget
if.target<-sapply(mat_pct_sigfig_manuscript_spec, function(x){
  x$tar
});names(if.target)<-names(mat_pct_sigfig_manuscript_spec)
# saveRDS(if.target, paste0(dir.o, "ifTarget_68.rds"))

### Plot
dir.o<-"/data_share/"
toplot<-readRDS(paste0(dir.o, "/signature.score_norm.med.rds"))
colnames(toplot)[grep("Smooth muscle cell",colnames(toplot))]<-"Smooth muscle (eMSC enriched)*"

# permutation results
permute_results<-readRDS(paste0(dir.o, "permute.verdict.rds"))

# ratioNext
spec_toplot<-readRDS(paste0(dir.o, "specificity.score_ratio.rds"))
rowpal_spec<-c("grey90", "#156D9F")
rowpal_spec<-colorRampPalette(rowpal_spec, space="Lab")(9)
rowcol_spec<-rowpal_spec[1+8*spec_toplot/max(spec_toplot)]
names(rowcol_spec)<-names(spec_toplot)

# onTarget
tar_toplot<-readRDS(paste0(dir.o, "ifTarget_68.rds"))
rowpal_tar<-c("gray60","grey30","#56B4E9")
rowcol_tar<-rowpal_tar[1+tar_toplot]
names(rowcol_tar)<-names(tar_toplot)

# permute test
palRow<-c("gray","black")

# row and column orders
library("cba")
coldist=as.dist(1-abs(cor(toplot,method="pearson")))
rowdist=as.dist(1-abs(cor(t(toplot),method="pearson")))
hc <- hclust(coldist, method="ward.D2")
hr <- hclust(rowdist, method="ward.D2")

optimal.row <- order.optimal(rowdist,hr$merge)
optimal.col <- order.optimal(coldist,hc$merge)
hr$merge <- optimal.row$merge; hr$order <- optimal.row$order
hc$merge <- optimal.col$merge; hc$order <- optimal.col$order

color.palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")),space="Lab")
palette.breaks<-seq(min(toplot),max(toplot),1)

library(gplots)
pdf(paste0(dir.plot,"pct_significant.genes_fc1_padj.05_manuscript_norm.med_publication.pdf"),width=10.2,height=3.4)
heatmap.2(as.matrix(t(toplot)), 
          col=color.palette,
          trace="none", scale="none",density.info="none",key.title="",
          breaks = palette.breaks,
          key.xlab="Median normed\nsignature score",
          #ColSideColors = rowcol_spec[rownames(toplot)], #ratioNext
          ColSideColors = rowcol_tar[rownames(toplot)], #onTarget
          key.par=list(mgp=c(1.5, 0.5, 0),
                       mar=c(5,1,5.5,1)),
          # key.xtickfun=function(){
          #   list(at=seq(0,1,0.25),labels<-seq(0,100,25))
          # },
          adjCol=1,cexCol = 0.9,cexRow = 1.1,srtCol=45,
          margins = c(9,16),lwid=c(0.2,6),lhei=c(1,4.5,4),
          dendrogram = "both",Colv = as.dendrogram(hr),Rowv = as.dendrogram(hc),
          colCol = palRow[permute_results[rownames(toplot)]+1]
)
dev.off()