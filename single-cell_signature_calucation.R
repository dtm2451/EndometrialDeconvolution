#* Author: Wanxin Wang
#* This script calculates the signature scores based on methods by
#* F. Pont, M. Tosolini, J. J. Fournié, Nucleic Acids Res. 47, e133 (2019).
#* for violinplots in Fig.3C, 5C,D, S7

rm(list=ls())
library(dplyr)
library(Seurat)
library(RColorBrewer)
library(plyr)
library(ggplot2)

dir.plot<-""
dir.o<-"/datashare/"

xcell_sc<-readRDS("GSE111976_ct_endo_10x.rds") # can be downloaded from GSE111976
xcell_sig_all_list<-readRDS(paste0(dir.o, "xcell_sig_all_list.rds"))
xcell_ctp<-names(xcell_sig_all_list)

# calculate signature score as the ratio between transcripts (UMI) that encode genes in the xCell signature to all transcripts (UMI) detected in each single cell (Methods)
tpm_total<-xcell_sc$nCount_RNA
tpm_total_check<-Matrix::colSums(xcell_sc@assays[["RNA"]]@counts)
plot(tpm_total,tpm_total_check)
rm("tpm_total_check")
ct<-xcell_sc@assays[["RNA"]]@counts
scores_by.sig<-lapply(1:length(xcell_sig_all_list), function(x){
  genes<-xcell_sig_all_list[[x]];print(paste0("l_genes: ",length(genes)))
  ct_sig<-ct[rownames(ct)%in%genes,];print(paste0("dim_genes: ",dim(ct_sig)))
  score<-Matrix::colSums(ct_sig)/tpm_total
});names(scores_by.sig)<-names(xcell_sig_all_list)

df_scores_by.sig<-as.data.frame(t(matrix(unlist(scores_by.sig),nrow=length(scores_by.sig), byrow=T)))
colnames(df_scores_by.sig)<-names(scores_by.sig)
rownames(df_scores_by.sig)<-names(scores_by.sig[[1]])
#saveRDS(df_scores_by.sig, paste0(dir.o, "df_scores_by.sig.rds"))


