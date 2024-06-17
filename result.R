library(dplyr)
library(tinyarray)
library(rtracklayer)
library(Seurat)

gtf <- rtracklayer::import('gencode.v22.annotation.gtf')
gtf <- as.data.frame(gtf)
dim(gtf)
head(gtf)
anno=dplyr::select(gtf,c(gene_name,gene_id))
anno<-distinct(anno,gene_id,.keep_all = T)


TCGA<-read.table("TCGA-THCA.htseq_counts.tsv",header = T,row.names = 1,sep = "\t",quote = "",check.names = F)#521
TCGA<-sam_filter(TCGA)#tinyarray包去除重复肿瘤样本,剩余560
t_index <- which(substr(colnames(TCGA),14,15) == '01')
tumor<-TCGA[,t_index]#456
normal<-TCGA[,-t_index]#41

#colnames(tumor) <- substr(colnames(tumor),1,12) #整理样本名，这里都为肿瘤，所以只取前12个字符
#colnames(normal) <- substr(colnames(normal),1,12) #整理样本名，这里都为肿瘤，所以只取前12个字符
TCGA<-cbind(tumor,normal)#560
#ID转换
TCGA$gene_id<-rownames(TCGA)
#TCGA$gene_id=unlist(str_split(TCGA$gene_id,"[.]",simplify=T))[,1]#去除小数点后的版本号
TCGA=merge(anno,TCGA,by ="gene_id")
TCGA<-na.omit(TCGA)
TCGA$gene_id<-NULL
TCGA<-TCGA[!duplicated(TCGA$gene_name),]
rownames(TCGA)<-TCGA$gene_name
TCGA$gene_name<-NULL
max(TCGA)
load("TE.Rdata")
cnv<-read.table("kmeans_df_s5.txt",header = T,sep = "\t")
cnv$infercnv[cnv$kmeans_class%in%c(2)]<-"para-tumors"
cnv$infercnv[!(cnv$kmeans_class%in%c(2))]<-"tumors"

scRNAsub.T@meta.data$inferCNV<-cnv[match(colnames(scRNAsub.T),rownames(cnv)),"infercnv"]
Idents(scRNAsub.T)<-"inferCNV"
markers <- FindAllMarkers(object = scRNAsub.T, only.pos = TRUE, min.pct = 0.25)
top = 50 # can be adjusted as needed
TopMarkers <- markers %>% filter(p_val_adj < 0.01) %>% group_by(cluster) %>%
  top_n(n = top, wt = avg_log2FC)
immunity<-TopMarkers[,c("cluster","gene")]
immunity<-immunity%>%split(.,.$cluster)%>%
  lapply(., function(x)(x$gene))
#每种细胞类型特征基因去重
immunity<-lapply(immunity, unique)
immunity[1]
library(GSVA)
TCGA_gsva<-as.data.frame(t(gsva(as.matrix(TCGA),immunity,method="ssgsea")))
#TCGA_gsva<-as.data.frame(t(TCGA_gsva))

library(ggplot2)
library(reshape2)
library(ggpubr)
library(tidyverse)
boxplot<-TCGA_gsva
boxplot$Type<-c(rep("Tumor",502),rep("Normal",58))
#boxplot<-rename(boxplot,Group=group)
melt <- melt(boxplot,
             id.vars = c("Type"),
             variable.name ="Signature",
             value.name = "Signature_score")
melt$Signature<-factor(melt$Signature,levels = c("para-tumors","tumors"))
box<- ggplot(melt,
             aes(x=Signature, y=Signature_score, 
                 fill = Type #??????????ɫ
             )) + #???????߿???ɫ
  geom_boxplot(notch = F, alpha = 0.95, 
               outlier.shape = 16,
               outlier.colour = "black", #outlier???ú?ɫ
               outlier.size = 0.65) +
  #?Զ?????ɫ
  scale_fill_manual(values= c("#377EB8","#E41A1C")) +
  ggtitle("K mean clustering coefficient 5") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5,size = 20,face = "bold"), 
        axis.text.y = element_text(angle = 0, size = 15,face = "bold"),
        axis.title.y = element_text(angle = 90, size = 20,face = "bold"),
        axis.title.x = element_text(angle = 0, size = 20,face = "bold"),
        legend.title = element_text(size = 20,face =  "bold"),
        title = element_text(size = 20,face =  "bold",hjust = 0.5),
        legend.text = element_text(size = 20,face = "bold")
  ) +
  xlab("")+ylab("Enrichment Score")+
  theme(legend.position = "top")

# ??עp value
box + stat_compare_means(label = "p.format",size=6)
ggsave("boxplot_5.pdf",width = 8,height =6)
###根据箱式图P值确定聚类系数K=11
load("TE.Rdata")
cnv<-read.table("kmeans_df_s11.txt",header = T,sep = "\t")
a<-scRNAsub.T@meta.data
cnv$group<-a[match(rownames(cnv),rownames(a)),"group"]
cnv$infercnv[cnv$kmeans_class%in%c(5,6)]<-"para-tumors"
cnv$infercnv[!(cnv$kmeans_class%in%c(5,6)) & cnv$group%in%"primary-tumors"]<-"primary-tumors"
cnv$infercnv[!(cnv$kmeans_class%in%c(5,6)) & cnv$group%in%"distant-metastases"]<-"distant-metastases"
cnv$infercnv[!(cnv$kmeans_class%in%c(5,6)) & cnv$group%in%"Lymph-Node"]<-"Lymph-Node"
cnv$infercnv[!(cnv$kmeans_class%in%c(5,6)) & cnv$group%in%"para-tumors"]<-"para-tumors-tumors"


scRNAsub.T@meta.data$inferCNV<-cnv[match(colnames(scRNAsub.T),rownames(cnv)),"infercnv"]


save(scRNAsub.T,file = "TE_inferCNV.Rdata")
