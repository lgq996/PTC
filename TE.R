#https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/

#安装发行版，作者推荐
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("infercnv")
#安装github上的最新版
#library("devtools")
#devtools::install_github("broadinstitute/infercnv")
BiocManager::install("infercnv",force = TRUE)
BiocManager::install("coda",force = TRUE)
BiocManager::install("rjags",force = TRUE)
library(coda)
library(rjags)
library(infercnv)
library(tidyverse)
library(Seurat)
library(SeuratObject)
load("TE.Rdata")
scRNAsub.T@meta.data$CNV<-scRNAsub.T@meta.data$group
scRNAsub.T@meta.data$CNV[scRNAsub.T@meta.data$CNV=="para-tumors"]<-"para-tumors(non-malignant)"
ann<-read.table("hg38_gencode_v27.txt",header = F,sep = "\t",quote = "")
##读取示例数据目录
exprMatrix =as.matrix(GetAssayData(scRNAsub.T,slot = "count")) 
cellAnnota = subset(scRNAsub.T@meta.data,select = "CNV")
geneLocate =ann[which(ann$V1%in%row.names(scRNAsub.T)),]
write.table(exprMatrix,"exprMatrix.txt",col.names = T,row.names = T,sep = "\t",quote = F)
write.table(cellAnnota,"cellAnnota.txt",col.names = F,row.names =T,sep = "\t",quote = F)
write.table(geneLocate,"geneLocate.txt",col.names = F,row.names = F,sep = "\t",quote = F)

#创建inferCNV对象，直接给相应的文件路径即可
infercnv_obj = CreateInfercnvObject(delim = '\t',
                                    raw_counts_matrix = "exprMatrix.txt",
                                    annotations_file = "cellAnnota.txt",
                                    gene_order_file = "geneLocate.txt",
                                    ref_group_names = c("para-tumors(non-malignant)"))

##分析细胞CNV
#cutoff阈值，Smart-seq2数据选“1”, 10x Genomics数据选“0.1”
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, 
                             out_dir='inferCNV', 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE)
 