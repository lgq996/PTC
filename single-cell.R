library(Seurat)
library(tidyverse)
library(patchwork)
library(ggsci)
mypal<-pal_npg("nrc")(10)

set.seed(123)  #设置随机数种子，使结果可重复
##==合并数据集==##
##使用目录向量合并
dir = c('data/PTC1P', 
        'data/PTC1T',
        'data/PTC2LeftLN',
        'data/PTC2P',
        'data/PTC2T',
        'data/PTC3LeftLN',
        'data/PTC3P',
        'data/PTC3RightLN',
        'data/PTC3T',
        'data/PTC4SC',
        'data/PTC5P',
        'data/PTC5RightLN',
        'data/PTC5T',
        'data/PTC6RightLN',
        'data/PTC7RightLN',
        'data/PTC8P',
        'data/PTC8T',
        'data/PTC9P',
        'data/PTC9T',
        'data/PTC10RightLN',
        'data/PTC10T',
        'data/PTC11RightLN',
        'data/PTC11SC')
names(dir) = c('PTC1P', 'PTC1T', 'PTC2LeftLN', 'PTC2P', 'PTC2T', 
               'PTC3LeftLN', 'PTC3P', 'PTC3RightLN', 'PTC3T', 'PTC4SC',
               'PTC5P', 'PTC5RightLN', 'PTC5T', 'PTC6RightLN', 'PTC7RightLN',
               'PTC8P', 'PTC8T', 'PTC9P', 'PTC9T', 'PTC10RightLN',
               'PTC10T', 'PTC11RightLN', 'PTC11SC')
#使用merge函数合并seurat对象
scRNAlist <- list()
#以下代码会把每个样本的数据创建一个seurat对象，并存放到列表scRNAlist里
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, min.cells = 3, min.features = 200)
}
#数据整合之前要对每个样本的seurat对象进行数据标准化和选择高变基因
# normalize and identify variable features for each dataset independently
scRNA <- lapply(X = scRNAlist, FUN = function(x) {
  x <- NormalizeData(x,normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = scRNA)
immune.anchors <- FindIntegrationAnchors(object.list = scRNA, anchor.features = features)
# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "RNA"

immune.combined[["percent.mt"]] <- PercentageFeatureSet(object = immune.combined, pattern = "^MT-")
pdf(file="04.featureViolin.pdf",width=10,height=6)          
VlnPlot(object = immune.combined, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)
dev.off()
immune.combined <- subset(x = immune.combined, 
                          subset = nFeature_RNA > 500 & nFeature_RNA<3000 & percent.mt < 10 & nCount_RNA>500 & nCount_RNA<5000)    


#TPM不用这行
immune.combined <- NormalizeData(object = immune.combined, normalization.method = "LogNormalize", scale.factor = 10000)
#修改方法
immune.combined <- FindVariableFeatures(object = immune.combined, selection.method = "vst", nfeatures = 5000)
#????????????ͼ
top10 <- head(x = VariableFeatures(object = immune.combined), 10)
pdf(file="04.featureVar.pdf",width=10,height=6)              #????????????????ͼ
plot1 <- VariableFeaturePlot(object = immune.combined,cols = mypal[1:2])
plot2 <- LabelPoints(plot = plot1, points = top10, xnudge = 0,ynudge = 0,
                     repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()

###################################05.PCA???ɷַ???###################################
##PCA????
immune.combined=ScaleData(immune.combined)                    
immune.combined=RunPCA(object= immune.combined,npcs = 20,pc.genes=VariableFeatures(object = immune.combined))    

#????ÿ??PCA?ɷֵ????ػ???
pdf(file="05.pcaGene.pdf",width=18,height=6)
VizDimLoadings(object = immune.combined, dims = 1:10, 
               reduction = "pca",nfeatures = 10,ncol = 5,col = mypal)
dev.off()

#???ɷַ???ͼ??
pdf(file="05.PCA.pdf",width=6.5,height=6)
DimPlot(object = immune.combined, reduction = "pca",cols = mypal[1])
dev.off()

#???ɷַ?????ͼ
pdf(file="05.pcaHeatmap.pdf",width=10,height=8)
DimHeatmap(object = immune.combined, dims = 1:10, cells = 500, fast = F,
           balanced = TRUE,nfeatures = 30,ncol=5)

dev.off()

#ÿ??PC??pֵ?ֲ??;??ȷֲ?
immune.combined <- JackStraw(object = immune.combined, num.replicate = 100)
immune.combined <- ScoreJackStraw(object = immune.combined, dims = 1:20)
pdf(file="05.pcaJackStraw.pdf",width=8,height=6)
JackStrawPlot(object = immune.combined, dims = 1:20,
              cols = colorRampPalette(colors = c(mypal))(20))
dev.off()





###################################06.TSNE??????????marker????###################################
##TSNE????????
pcSelect=20
immune.combined <- FindNeighbors(object = immune.combined, dims = 1:pcSelect)
#先粗略聚类
immune.combined <- FindClusters(object = immune.combined, resolution = 0.1)###粗聚类0.1                 
immune.combined <- RunTSNE(object = immune.combined, dims = 1:pcSelect)                      
pdf(file="06.TSNE.pdf",width=8,height=6)
DimPlot(object = immune.combined, label = TRUE,pt.size = 1,label.size = 10,label.box = T,repel =T,group.by = "seurat_clusters",
        cols = colorRampPalette(colors = c(mypal))(14))+NoLegend()+
  labs(x="tSNE1",y="tSNE2",title = "")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 15,face = "bold"), 
        axis.text.y = element_text(angle = 0, size = 15,face = "bold"),
        axis.title.y = element_text(angle = 90, size = 20,face = "bold"),
        axis.title.x = element_text(angle = 0, size = 20,face = "bold"),
        legend.title = element_text(size = 20,face =  "bold"),
        legend.text = element_text(size = 20,face = "bold"),
        panel.border = element_rect(fill = NA,colour = "black",size = 2,linetype = "solid"))
dev.off()
write.table(immune.combined$seurat_clusters,file="06.tsneCluster_0.1.txt",quote=F,sep="\t",col.names=F)
###根据marker基因进行细胞类型注释
features<-c("LYZ", "S100A8", "S100A9", "CD14",
            "CD3D", "CD3E","IL7R","IL32","TRAC",
            "CD79A", "CD79B","MS4A1","CD74",
            "TG","CLU","FN1","MGST1","S100A13",
            "RGS5","IGFBP7","TAGLN","COL1A2","ACTA2",
            "TIMP3","RAMP2","CLDN5","TFPI","MGP")
class<-c("Myeloid cells","Myeloid cells","Myeloid cells","Myeloid cells",
         "T/NK cells","T/NK cells","T/NK cells","T/NK cells","T/NK cells",
         "B cells","B cells","B cells","B cells",
         "Thyroid epithelial cells","Thyroid epithelial cells","Thyroid epithelial cells","Thyroid epithelial cells","Thyroid epithelial cells",
         "Fibroblasts","Fibroblasts","Fibroblasts","Fibroblasts","Fibroblasts",
         "Endothelial cells","Endothelial cells","Endothelial cells","Endothelial cells","Endothelial cells")
p1<-DotPlot(immune.combined,features = features,cols = mypal,group.by = "cellType_marker")
pdot<-p1$data
pdot$class<-rep(class,6)
pdot$features.plot<-factor(pdot$features.plot,levels = unique(pdot$features.plot))
pdot$class<-factor(pdot$class,levels = c("Myeloid cells","T/NK cells","B cells",
                                                  "Thyroid epithelial cells","Fibroblasts","Endothelial cells"))
pdot$id<-factor(pdot$id,levels = c("Myeloid cells","T/NK cells","B cells",
                                         "Thyroid epithelial cells","Fibroblasts","Endothelial cells"))

pdot1<-ggplot(pdot,aes(x=features.plot,y=id))+
  geom_point(aes(fill=avg.exp,size=pct.exp),color='black',shape=21)+
  theme_bw(base_size = 14)+
  xlab("")+ylab("")+
  scale_fill_gradient2(low = '#4DBBD5FF',mid = 'white',high = '#E64B35FF',midpoint = 2,name='Mean expression')+
  scale_size(range = c(1,6),name='Percentage')+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        aspect.ratio = 0.5,
        plot.margin = margin(t=1,r=1,b=1,l=1,unit = 'cm'),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,face = 'bold'),
        axis.text.y = element_text(angle = 0,hjust = 1,vjust = 0.5,face = 'bold'))+
  coord_cartesian(clip = 'off')
pdot_test<-pdot1+
  scale_x_discrete(position = 'top')+
  theme(axis.text.x = element_text(angle = 90,hjust = 0))
library(jjAnno)
annoRect(object = pdot_test,annoPos = "top",aesGroup = T,aesGroName = 'class',yPosition = c(6.6,9.2),
         rectWidth = 1,
         rectAngle = 90,normRectShift = 0.1,textHVjust =1.5,textRot = 45,hjust = 0,textShift = 0,
         rotateRect = T,alpha = 0.5,
         addText = T)
ggsave("cell_bubble_celltype.pdf",width =10,height =8)
#确定好后在做决定
immune.combined@meta.data$cellType_marker<-ifelse(immune.combined@meta.data$seurat_clusters%in%c(0,1,4),"T/NK cells",
                                       ifelse(immune.combined@meta.data$seurat_clusters%in%c(3,8),"B cells",
                                              ifelse(immune.combined@meta.data$seurat_clusters%in%c(2,5,9,10,13),"Thyroid epithelial cells",
                                                     ifelse(immune.combined@meta.data$seurat_clusters%in%c(7),"Fibroblasts",
                                                            ifelse(immune.combined@meta.data$seurat_clusters%in%c(6),"Myeloid cells","Endothelial cells")))))
pdf(file="06.TSNE_cellType_marker.pdf",width=8,height=6)
DimPlot(object = immune.combined, label = TRUE,pt.size = 1,label.size = 6,label.box = T,repel =T,group.by = "cellType_marker",
        cols = mypal[1:6])+NoLegend()+
  labs(x="tSNE1",y="tSNE2",title = "")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 15,face = "bold"), 
        axis.text.y = element_text(angle = 0, size = 15,face = "bold"),
        axis.title.y = element_text(angle = 90, size = 20,face = "bold"),
        axis.title.x = element_text(angle = 0, size = 20,face = "bold"),
        legend.title = element_text(size = 20,face =  "bold"),
        legend.text = element_text(size = 20,face = "bold"),
        panel.border = element_rect(fill = NA,colour = "black",size = 2,linetype = "solid"))
dev.off()

immune.combined@meta.data$HT<-ifelse(immune.combined@meta.data$orig.ident%in%c("PTC1P","PTC1T",
                                                                               "PTC2LeftLN","PTC2P","PTC2T",
                                                                               "PTC4SC",
                                                                               "PTC8P","PTC8T"),"HT","Non-HT")
immune.combined@meta.data$group<-ifelse(immune.combined@meta.data$orig.ident%in%c("PTC1P","PTC2P","PTC3P","PTC5P","PTC8P","PTC9P"),"para-tumors",
                                        ifelse(immune.combined@meta.data$orig.ident%in%c("PTC1T","PTC2T","PTC3T","PTC5T","PTC8T","PTC9T","PTC10T"),"primary-tumors",
                                               ifelse(immune.combined@meta.data$orig.ident%in%c("PTC2LeftLN","PTC3LeftLN","PTC3RightLN","PTC5RightLN","PTC6RightLN","PTC7RightLN","PTC10RightLN","PTC11RightLN"),"Lymph-Node","distant-metastases")))
mypal1<-pal_lancet("lanonc")(9)

###展示各个分类的降维图
library(egg)##控制绘图区大小，以保证标签文字等所占空间不同时，同批图像还是同样大小
library(grid)
aaa<-DimPlot(object = immune.combined, pt.size = 1,label=F,
             label.size = 10,label.box = T,repel =T,group.by = c("orig.ident"),
             cols = colorRampPalette(colors = c(mypal1[1:8]))(23))+
  labs(x="tSNE1",y="tSNE2",title = "",colour="orig.ident")+
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title.y = element_text(angle = 90, size = 20,face = "bold"),
        axis.title.x = element_text(angle = 0, size = 20,face = "bold"),
        legend.title = element_text(size = 20,face =  "bold"),
        legend.text = element_text(size = 15,face = "bold"),
        panel.border = element_blank())
ggsave(file="06.TSNE_orig.ident.pdf",egg::set_panel_size(aaa,width=unit(4.5, "in"), height=unit(5, "in")), 
       width =12, height = 6, units = 'in', dpi = 300)

aaa<-DimPlot(object = immune.combined, pt.size = 1,label=F,
             label.size = 10,label.box = T,repel =T,group.by = c("group"),
             cols = mypal1[1:4])+
  labs(x="tSNE1",y="tSNE2",title = "",colour="group")+
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title.y = element_text(angle = 90, size = 20,face = "bold"),
        axis.title.x = element_text(angle = 0, size = 20,face = "bold"),
        legend.title = element_text(size = 20,face =  "bold"),
        legend.text = element_text(size = 15,face = "bold"),
        panel.border = element_blank())
ggsave(file="06.TSNE_group.pdf",egg::set_panel_size(aaa,width=unit(4.5, "in"), height=unit(5, "in")), 
       width = 8, height = 6, units = 'in', dpi = 300)

aaa<-DimPlot(object = immune.combined, pt.size = 1,label=F,
             label.size = 10,label.box = T,repel =T,group.by = c("HT"),
             cols = mypal1[4:5])+
  labs(x="tSNE1",y="tSNE2",title = "",colour="HT")+
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title.y = element_text(angle = 90, size = 20,face = "bold"),
        axis.title.x = element_text(angle = 0, size = 20,face = "bold"),
        legend.title = element_text(size = 20,face =  "bold"),
        legend.text = element_text(size = 15,face = "bold"),
        panel.border = element_blank())
ggsave(file="06.TSNE_HT.pdf",egg::set_panel_size(aaa,width=unit(4.5, "in"), height=unit(5, "in")), 
       width = 8, height = 6, units = 'in', dpi = 300)

###百分比堆积条形图
bar<-data.frame(immune.combined@meta.data)
bar<-bar[,c(7:9)]
bar<-reshape2::dcast(bar,HT+group~cellType_marker,fill = 0)
bar<-reshape2::melt(bar)
#bar$group<-factor(bar$group,levels = c("HT","Non-HT"))
ggplot(data = bar,aes(x=HT,y=value))+
  geom_bar(aes(fill=variable),stat = 'identity',position = position_fill(),width = 0.5)+
  theme_bw()+scale_fill_manual(values =c("B cells"=mypal[1],"Endothelial cells"=mypal[2],
                                         "Fibroblasts"=mypal[3],"T/NK cells"=mypal[4],
                                         "Thyroid epithelial cells"=mypal[5],
                                         "Myeloid cells"=mypal[6]))+
  facet_grid(~group,scales = "free")+
  labs(x="",y="",title = "")+labs(fill="CellType")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 15,face = "bold"), 
        axis.text.y = element_text(angle = 0, size = 20,face = "bold"),
        axis.title.y = element_text(angle = 90, size = 20,face = "bold"),
        axis.title.x = element_text(angle = 0, size = 20,face = "bold"),
        legend.title = element_text(size = 20,face =  "bold"),
        legend.text = element_text(size = 8,face = "bold"),
        panel.border = element_blank(),
        strip.text.x = element_text(angle = 0, size = 10,face = "bold"),
        strip.text.y = element_text(angle = 270, size = 15,face = "bold"),
        strip.background.x = element_rect(
          color="white", fill=mypal[5], size=1.5, linetype="solid"
        ),
        strip.background.y = element_rect(
          color="white", fill=mypal[2], size=1.5, linetype="solid"
        ))
ggsave("bar_all.pdf",width =10,height = 6)
###用药前后各种免疫细胞的差异条形图
#思路是先不分区域来整体比较每种细胞类型的差异是否显著，选择差异显著的细胞类型进行进一步分析
library(rstatix)
library(SeuratDisk)
cellinfo <- Seurat::FetchData(immune.combined, vars = c("HT", "group","cellType_marker"))
#cellinfo<-cellinfo[which(cellinfo$group=="distant-metastases"),]
cell.col <- setNames(object = c(mypal[1], mypal[2]),
                     nm = c("HT", "Non-HT"))

# 构建seurat_clusters和Group的列联表
tbl <- table(cellinfo$cellType_marker, cellinfo$HT)
tbl
# 进行fisher精确检验和post hoc检验
# post hoc检验：对各组样本进行 one vs other的fisher检验，进行多重性校正，得到各组的p-adj
fisher.test(tbl, simulate.p.value = T) # fisher精确检验
post.hoc <- row_wise_fisher_test(tbl) # post hoc检验
post.hoc$sig.label <- ifelse(test = post.hoc$p.adj.signif == "ns", # 调整显著性显示标签
                             yes = "", no = post.hoc$p.adj.signif) # 不显示NS (No Significant)
str(post.hoc)
head(post.hoc)
# 绘制柱状图
plot.data <- as.data.frame(tbl)
head(plot.data)
p2 <- ggplot(plot.data, aes(x = Var1, y = Freq, fill = Var2)) + 
  geom_bar(stat = "identity", position = position_fill(reverse = T)) +  # 绘制堆叠比例柱状图，将WT和KO的顺序倒过来
  scale_fill_manual(values = cell.col) +                                # 设置不同组对应的颜色
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +                   # 设置y坐标轴的刻度
  labs(title = "",fill="")+
  geom_text(data = post.hoc,                                            # 设置顶部的显著性标签
            aes(x = group, y = 1.1, label = sig.label),                 # 标签的位置和具体内容
            inherit.aes = F) +                                          
  xlab("") +                                                     # 修改x和y轴列名
  ylab("Fraction of Cells") +
  theme_classic() +                                                     # 去除不必要的元素
  theme(plot.title =element_text(size = 25,face = "bold",hjust = 0.5,colour = "red"),                           
        axis.text = element_text(size = 10,face = "bold"),              # 调整坐标轴刻度字号大小
        axis.text.x = element_text(face = "bold",angle = 45,hjust = 1),
        axis.title = element_text(size = 20,face = "bold"),
        legend.title = element_text(size = 20,face =  "bold"),          # 调整坐标轴标题字号大小
        legend.text = element_text(size = 10,face = "bold"),
  )                           
p2
ggsave(filename = "ALL.pdf", width = 8, height = 6)
###分别提取每个细胞亚群
table(immune.combined@meta.data$cellType_marker)
Cells.sub.T <- subset(immune.combined@meta.data, cellType_marker=="T/NK cells")
scRNAsub.T <- subset(immune.combined, cells =row.names(Cells.sub.T))
save(scRNAsub.T,file = "T.Rdata")

Cells.sub.T <- subset(immune.combined@meta.data, cellType_marker=="B cells")
scRNAsub.T <- subset(immune.combined, cells =row.names(Cells.sub.T))
save(scRNAsub.T,file = "B.Rdata")

Cells.sub.T <- subset(immune.combined@meta.data, cellType_marker=="Endothelial cells")
scRNAsub.T <- subset(immune.combined, cells =row.names(Cells.sub.T))
save(scRNAsub.T,file = "Endothelial.Rdata")

Cells.sub.T <- subset(immune.combined@meta.data, cellType_marker=="Fibroblasts")
scRNAsub.T <- subset(immune.combined, cells =row.names(Cells.sub.T))
save(scRNAsub.T,file = "Fibroblasts.Rdata")

Cells.sub.T <- subset(immune.combined@meta.data, cellType_marker=="Thyroid epithelial cells")
scRNAsub.T <- subset(immune.combined, cells =row.names(Cells.sub.T))
save(scRNAsub.T,file = "TE.Rdata")

Cells.sub.T <- subset(immune.combined@meta.data, cellType_marker=="Myeloid cells")
scRNAsub.T <- subset(immune.combined, cells =row.names(Cells.sub.T))
save(scRNAsub.T,file = "Myeloid.Rdata")

#以下没有保存

###细胞通讯
#通讯分析
library(Seurat)
library(tidyverse)
library(iTALK)
library(dplyr)
a<-read.table("CD4_Tn_Treg.txt",header = F,sep = "\t",quote = "")
b<-read.table("CD8_Tex_Teff.txt",header = F,sep = "\t",quote = "")
c<-read.table("cDC.txt",header = F,sep = "\t",quote = "")
d<-read.table("Endothelial_13.txt",header = F,sep = "\t",quote = "")
e<-read.table("Fibroblasts_4.txt",header = F,sep = "\t",quote = "")
f<-read.table("Intermediate_B.txt",header = F,sep = "\t",quote = "")
g<-read.table("mTE3.txt",header = F,sep = "\t",quote = "")
h<-read.table("nTE0.txt",header = F,sep = "\t",quote = "")
i<-read.table("nTE2.txt",header = F,sep = "\t",quote = "")
sub<-subset(immune.combined, cells =c(a$V1,b$V1,c$V1,d$V1,e$V1,f$V1,g$V1,h$V1,i$V1))
sub@meta.data$iTALK<-ifelse(colnames(sub)%in%a$V1,"CD4_Tn_Treg",
                  ifelse(colnames(sub)%in%b$V1,"CD8_Tex_Teff",
                         ifelse(colnames(sub)%in%c$V1,"cDC",
                                ifelse(colnames(sub)%in%d$V1,"Endothelial_13",
                                       ifelse(colnames(sub)%in%e$V1,"Fibroblasts_4",
                                              ifelse(colnames(sub)%in%f$V1,"Intermediate_B",
                                                     ifelse(colnames(sub)%in%g$V1,"mTE3",
                                                            ifelse(colnames(sub)%in%h$V1,"nTE0","nTE2"))))))))
##准备iTALK输入文件
# 输入文件格式为数据框，行为细胞名称，列为基因名和细胞注释信息，基因名命名的列值为表达数据。
# 细胞注释信息必选“cell_type”，值为细胞类型注释信息；可选“compare_group”，值为细胞的样本分组信息。
###HT样本的T(患者1,2,8)
#ht1<-subset(x=sub,orig.ident%in%c("PTC1P","PTC1T"))#1和2的细胞数太少，此处用8的
#ht2<-subset(x=sub,orig.ident%in%c("PTC2LeftLN","PTC2P","PTC2T"))
#ht8<-subset(x=sub,orig.ident%in%c("PTC8P","PTC8T"))
###非HT样本的T(患者3,5,9)
#nht3<-subset(x=sub,orig.ident%in%c("PTC3LeftLN","PTC3P","PTC3RightLN","PTC3T"))
#nht5<-subset(x=sub,orig.ident%in%c("PTC5P","PTC5RightLN","PTC5T"))
#nht9<-subset(x=sub,orig.ident%in%c("PTC9P","PTC9T"))

##==细胞通讯关系概览==##
cell.meta <- subset(sub@meta.data, select=c("orig.ident","iTALK"))
names(cell.meta) <- c("compare_group", "cell_type")
cell.expr <- data.frame(t(as.matrix(sub@assays$RNA@counts)), check.names=F)
data.italk <- merge(cell.expr, cell.meta, by=0)
rownames(data.italk) <- data.italk$Row.names
data.italk <- data.italk[,-1]
##设置绘图颜色和其他变量
mycolor <- c("#92D0E6","#F5949B","#E11E24","#FBB96F","#007AB8","#A2D184","#00A04E","#BC5627","#0080BD","#EBC379","#A74D9D")
mycolor<-jjAnno::useMyCol("stallion",n=20)
cell_type <- unique(data.italk$cell_type)
cell_col <- structure(mycolor[1:length(cell_type)], names=cell_type)
#通讯类型变量
comm_list<-c('growth factor','other','cytokine','checkpoint')
#圈图展示配体-受体对的数量
PN=20

paitent<-"All"
for (i in paitent) {
  dir.create(paste0('./',i))
  setwd(paste0('./',i))
  #实际上一个样本之内的细胞才能通讯，这里提取一组样本分析是为了克服单细胞数据稀疏性造成的误差
  #所有样本都应分析,找出共性的变化
  #data1 <- subset(data.italk, subset=data.italk$compare_group==i)#此处我们想看的就是肿瘤和癌旁，所以不分
  data1<-data.italk
  ##寻找高表达的配体受体基因,top_genes=50代表提取平均表达量前50%的基因
  highly_exprs_genes <- rawParse(data1, top_genes=50, stats="mean")
  res<-NULL
  for(comm_type in comm_list){
    #comm_type='other' 
    #多个细胞类型之间显著表达的配体-受体，结果会过滤同一细胞类型内的配体-受体关系
    res_cat <- FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
    #单个细胞类型显著表达的配体-受体
    #res_cat <- FindLR(data_1=deg_cd4T, datatype='DEG', comm_type=comm_type)
    res_cat <- res_cat[order(res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs, decreasing=T),]
    write.csv(res_cat, paste0('LRpairs_Overview_',comm_type,'.xls'))
    #plot by ligand category
    pdf(paste0('LRpairs_Overview_',comm_type,'_circ.pdf'), width=6, height=6)
    if(nrow(res_cat)==0){
      next
    }else if(nrow(res_cat)>=PN){
      LRPlot(res_cat[1:PN,], datatype='mean count', link.arr.lwd=res_cat$cell_from_mean_exprs[1:PN],
             cell_col=cell_col, link.arr.width=res_cat$cell_to_mean_exprs[1:PN])
    }else{
      LRPlot(res_cat, datatype='mean count', link.arr.lwd=res_cat$cell_from_mean_exprs,
             cell_col=cell_col, link.arr.width=res_cat$cell_to_mean_exprs)
    }
    dev.off()
    pdf(paste0('LRpairs_Overview_',comm_type,'_net.pdf'), width=6, height=6)
    NetView(res_cat, col=cell_col, vertex.label.cex=1.2, edge.label.cex=0.9, 
            vertex.size=30, arrow.width=3, edge.max.width=10, margin = 0.2)
    dev.off()
    res<-rbind(res,res_cat)
  }
  
  ##所有配体-受体分类一起作图
  res<-res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs, decreasing=T),]
  write.csv(res, 'LRpairs_Overview_all.xls')
  if(is.null(res)){
    print('No significant pairs found')
  }else if(nrow(res)>=PN){
    res.top <- res[1:PN,]
    pdf(paste0('LRpairs_Overview_all_net.pdf'), width=6, height=6)
    NetView(res, col=cell_col, vertex.label.cex=1.2, edge.label.cex=0.9, 
            vertex.size=30, arrow.width=3, edge.max.width=10, margin = 0.2)
    dev.off()    
    pdf('LRpairs_Overview_all_circ.pdf', width=6, height=6)
    LRPlot(res.top, datatype='mean count', link.arr.lwd=res.top$cell_from_mean_exprs,
           cell_col=cell_col, link.arr.width=res.top$cell_to_mean_exprs)
    dev.off()
  }else{
    pdf(paste0('LRpairs_Overview_all_net.pdf'), width=6, height=6)
    NetView(res, col=cell_col, vertex.label.cex=1.2, edge.label.cex=0.9, 
            vertex.size=30, arrow.width=3, edge.max.width=10, margin = 0.2)
    dev.off()    
    pdf('LRpairs_Overview_all_circ.pdf', width=6, height=6)
    LRPlot(res, datatype='mean count', link.arr.lwd=res$cell_from_mean_exprs,
           cell_col=cell_col, link.arr.width=res$cell_to_mean_exprs)
    dev.off()
  }  
  setwd("E:/lgq/重新分析")
}
###cellchat
devtools::install_github("sqjin/CellChat")
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(SeuratData)
options(stringsAsFactors = FALSE)
data.input  <- sub@assays$RNA@data
identity = data.frame(group =sub$iTALK   , row.names = names(sub$iTALK)) # create a dataframe consisting of the cell labels
unique(identity$group) # check the cell labels
cellchat <- createCellChat(object=data.input)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels")
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
CellChatDB <- CellChatDB.human 
(out3 <- capture.output(str(CellChatDB)))
out4 <- paste(out3, collapse="\n")
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
cellchat@DB <- CellChatDB.use # set the used database in the object
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
#future::plan("multiprocess", workers = 4) # do parallel  这里似乎有一些bug，在Linux上居然不行。de了它。
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human) 


cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
pdf("circle_all.pdf",width = 12,height = 6)
par(mfrow = c(1,2), xpd=T)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,   
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

#可以不用
mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
cellchat@netP$pathways
levels(cellchat@idents) 

netVisual_bubble(cellchat, sources.use = c(7:9), targets.use = c(1:6), remove.isolate = FALSE)
#根据上一步的热图确定了MIF通路为细胞间通讯的主导通路
pathways.show <- c("MIF") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(7,9) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout = "hierarchy")
#par(mfrow=c(1,1))
#netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

netAnalysis_contribution(cellchat, signaling = pathways.show)
##可视化单个配体-受体对
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(7,9) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  
                     pairLR.use = LR.show, vertex.receiver = vertex.receiver,layout = "hierarchy")
#使用小提琴/点图绘制信号基因表达分布
plotGeneExpression(cellchat, signaling = "MIF")
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 6, font.size = 10)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("MIF"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2
