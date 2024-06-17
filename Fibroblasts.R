load(file = "Fibroblasts.Rdata")
library(jjAnno)
library(Seurat)
library(Matrix)
library(stringr)
# CreateSeuratObject
library(Seurat)
library(limma)
#library(SingleR)
library(Seurat)
library(dplyr)
library(magrittr)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(ggsci)
library(scRNAtoolVis)
library(ggrepel)
library(reshape2)
library(BuenColors)
mypal<-pal_npg("nrc")(10)

#memory.limit(size = 100000)

scRNAsub.T <- FindVariableFeatures(object = scRNAsub.T, selection.method = "vst", nfeatures = 5000)

scRNAsub.T=ScaleData(scRNAsub.T)                    
scRNAsub.T=RunPCA(object= scRNAsub.T,npcs = 20,pc.genes=VariableFeatures(object = scRNAsub.T))    

scRNAsub.T <- JackStraw(object = scRNAsub.T, num.replicate = 100)
scRNAsub.T <- ScoreJackStraw(object = scRNAsub.T, dims = 1:20)
JackStrawPlot(object = scRNAsub.T, dims = 1:20,
              cols = colorRampPalette(colors = c(mypal))(20))
###################################06.TSNE??????????marker????###################################
##TSNE????????
pcSelect=20
scRNAsub.T <- FindNeighbors(object = scRNAsub.T, dims = 1:pcSelect)
#先粗略聚类
scRNAsub.T <- FindClusters(object = scRNAsub.T, resolution = 0.4)###粗聚类0.3，细聚类0.8                 
scRNAsub.T <- RunTSNE(object = scRNAsub.T, dims = 1:pcSelect)                      
pdf(file="06.TSNE_Fcell.pdf",width=8,height=6)
DimPlot(object = scRNAsub.T, label = TRUE,pt.size = 1,label.size = 10,label.box = T,repel =T,group.by = "seurat_clusters",
        cols = colorRampPalette(colors = c(mypal))(11))+NoLegend()+
  labs(x="tSNE1",y="tSNE2",title = "")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 15,face = "bold"), 
        axis.text.y = element_text(angle = 0, size = 15,face = "bold"),
        axis.title.y = element_text(angle = 90, size = 20,face = "bold"),
        axis.title.x = element_text(angle = 0, size = 20,face = "bold"),
        legend.title = element_text(size = 20,face =  "bold"),
        legend.text = element_text(size = 20,face = "bold"),
        panel.border = element_rect(fill = NA,colour = "black",size = 2,linetype = "solid"))
dev.off()
###根据marker基因进行细胞类型注释
#with IL-6, Ly6C and PDGFRα 
#marking inflammatory CAFs; Cxcl12 marking immune-regulatory CAFs; 
#MHC-II (H2-Aa, H2-Ab1, and CD74 in mice; HLA-DRA, HLA-DPA1, HLA-DQA, and CD74 inhuman) 
#marking antigen-presenting CAFs; and ACTA2, TAGLN and POSTN marking 
#myCAFs
features<-c("IL6","LY6G6C","PDGFRA",
            "CXCL12",
            "HLA-DRA","HLA-DPA1","HLA-DQA1","CD74",
            "ACTA2","TAGLN","POSTN")
class<-c("iCAFs","iCAFs","iCAFs",
         "irCAFs",
         "apCAFs","apCAFs","apCAFs","apCAFs",
         "myCAFs","myCAFs","myCAFs")
p1<-DotPlot(scRNAsub.T,features = features,cols = colorRampPalette(colors = c(mypal))(11),group.by = "seurat_clusters")
pdot<-p1$data
pdot$class<-rep(class,11)
pdot$features.plot<-factor(pdot$features.plot,levels = unique(pdot$features.plot))
pdot$class<-factor(pdot$class,levels = unique(pdot$class))

pdot<-ggplot(pdot,aes(x=features.plot,y=id))+
  geom_point(aes(fill=avg.exp,size=pct.exp),color='black',shape=21)+
  theme_bw(base_size = 10)+
  xlab("")+ylab("")+
  scale_fill_gradient2(low = '#4DBBD5FF',mid = 'white',high = '#E64B35FF',midpoint = 1,name='Mean expression')+
  scale_size(range = c(0.2,8),name='Percentage')+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        aspect.ratio = 0.5,
        #plot.margin = margin(t=1,r=1,b=1,l=1,unit = 'cm'),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,face = 'bold'),
        axis.text.y = element_text(angle = 0,hjust = 1,vjust = 0.5,face = 'bold'))+
  coord_cartesian(clip = 'off')
pdot_test<-pdot+
  scale_x_discrete(position = 'top')+
  theme(axis.text.x = element_text(angle = 90,hjust = 0))
annoRect(object = pdot_test,annoPos = "top",aesGroup = T,aesGroName = 'class',yPosition = c(11.6,15),rectWidth = 1,
         rectAngle = 90,normRectShift = 0,textHVjust =2.5,textRot = 45,hjust = 0,textShift = 0,
         rotateRect = T,alpha = 0.5,
         addText = T)
ggsave("cell_bubble_CAF.pdf",width = 8,height =8)
###通过富集分析进一步确定亚型
Idents(scRNAsub.T)<-"seurat_clusters"
markers.CAF.cluster<- FindAllMarkers(object = scRNAsub.T,
                                  only.pos = FALSE,
                                  min.pct = 0.25,
                                  logfc.threshold = 0)
top = 50 # can be adjusted as needed
c.TopMarkers <- markers.CAF.cluster %>% filter(p_val_adj < 0.01 & abs(avg_log2FC)>0.5) %>% group_by(cluster) %>%
  top_n(n = top, wt = avg_log2FC)
###根据每个cluster的marker基因可以做富集分析，来确定亚群的生物学作用
###绘制热图
library(ClusterGVis)
library(org.Hs.eg.db)
# prepare data from seurat object
st.data <- prepareDataFromscRNA(object = scRNAsub.T,
                                diffData = c.TopMarkers,
                                showAverage = F, keep.uniqGene = FALSE,
                                sep = "_")
library(R.utils)

R.utils::setOption("clusterProfiler.download.method",'auto')
# enrich for clusters
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 5,
                        seed = 5201314)
col<-c(rep(jjAnno::useMyCol("stallion",n = 11),each = 5))
# add GO annotation
pdf('GO_BP.pdf',height = 16,width = 20,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           column_title_rot = 90,
           markGenes = NULL,#基因数太多，此处不显示
           markGenes.side = "left",
           annoTerm.data = enrich,
           show_column_names = F,
           line.side = "left",
           show_row_dend = F,
           cluster.order = c(1:11),
           go.col = col,
           add.bar = T)
dev.off()





###此处不作细胞类型注释，而是要看每个cluster的分化过程
scRNAsub.T@meta.data$cellType_CAF<-ifelse(scRNAsub.T@meta.data$seurat_clusters%in%c(0,1,2,3,5,6),"myCAFs ",
                                          ifelse(scRNAsub.T@meta.data$seurat_clusters%in%c(4,8),"iCAFs_irCAFs_myCAFs",
                                                 ifelse(scRNAsub.T@meta.data$seurat_clusters%in%c(10),"irCAFs_apCAFs_myCAFs",
                                                        ifelse(scRNAsub.T@meta.data$seurat_clusters%in%c(9),"apCAFs_myCAFs","iCAFs_irCAFs_apCAFs_myCAFs"))))
scRNAsub.T$cellType_CAF<-NULL
#先画一个，获取图例
pdf(file="06.TSNE_CAFcell_cellType.pdf",width=8,height=6)
DimPlot(object = scRNAsub.T, label = T,pt.size = 1,label.size =6,label.box = T,repel =T,group.by = "cellType_CAF",
        cols = colorRampPalette(colors = c(mypal))(5))+NoLegend()+
  labs(x="tSNE1",y="tSNE2",title = "")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 15,face = "bold"), 
        axis.text.y = element_text(angle = 0, size = 15,face = "bold"),
        axis.title.y = element_text(angle = 90, size = 20,face = "bold"),
        axis.title.x = element_text(angle = 0, size = 20,face = "bold"),
        legend.title = element_text(size = 20,face =  "bold"),
        legend.text = element_text(size = 8,face = "bold"),
        panel.border = element_rect(fill = NA,colour = "black",size = 2,linetype = "solid"))
dev.off()

###画细胞比例堆积条形图
#devtools::install_github("caleblareau/BuenColors")
library(reshape2)
library(BuenColors)
jdb_palette("lawhoops")
bar<-scRNAsub.T@meta.data
bar<-bar[,c(8,9,6)]
bar<-reshape2::dcast(bar,HT+group~seurat_clusters,fill = 0)
bar<-melt(bar)
#bar$Group<-factor(bar$Group,levels = c("Pre","Post"))
ggplot(data = bar,aes(x=HT,y=value))+
  geom_bar(aes(fill=variable),stat = 'identity',position = position_fill(),width = 0.5)+
  theme_bw()+scale_fill_manual(values =c(colorRampPalette(colors = c(jdb_palette("lawhoops")))(11)))+
  #facet_grid(~group,scales = "free")+
  labs(x="",y="",title = "")+labs(fill="Cluster")+
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
ggsave("bar_CAF_all.pdf",width =6,height = 6)
library(rstatix)
library(SeuratDisk)
cellinfo <- Seurat::FetchData(scRNAsub.T, vars = c("HT", "group","seurat_clusters"))
#cellinfo<-cellinfo[which(cellinfo$group=="metastasis"),]
cell.col <- setNames(object = c(mypal[7], mypal[8]),
                     nm = c("HT", "Non-HT"))

# 构建seurat_clusters和Group的列联表
tbl <- table(cellinfo$seurat_clusters, cellinfo$HT)
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
ggsave(filename = "CAFbar_sig.pdf", width = 8, height = 6)
###HT样本的T(患者1,2,8)
ht1<-subset(x=scRNAsub.T,orig.ident%in%c("PTC1P","PTC1T"))#1和2的细胞数太少，此处用8的
ht2<-subset(x=scRNAsub.T,orig.ident%in%c("PTC2LeftLN","PTC2P","PTC2T"))
ht8<-subset(x=scRNAsub.T,orig.ident%in%c("PTC8P","PTC8T"))
###非HT样本的T(患者3,5,9)
nht3<-subset(x=scRNAsub.T,orig.ident%in%c("PTC3LeftLN","PTC3P","PTC3RightLN","PTC3T"))
nht5<-subset(x=scRNAsub.T,orig.ident%in%c("PTC5P","PTC5RightLN","PTC5T"))
nht9<-subset(x=scRNAsub.T,orig.ident%in%c("PTC9P","PTC9T"))
###轨迹分析
library(Seurat)
library(dplyr)
library(monocle)
library(ggsci)
#mypal<-pal_jco("default", alpha = 0.6)(9)
###此处要合并HT的1,2,8和nHT的3,5,9,与TE不同，这里不能分开分析，因为有些患者并不是都含有所有的细胞类型
###在分析基因随着拟时序变化时，要考虑特定细胞亚群的特征基因
pbmc <- merge(ht1, y = c(ht2,ht8), add.cell.ids = c("ht1", "ht2","ht8"), project = "HT",merge.data = TRUE)
pbmc <- merge(nht3, y = c(nht5,nht9), add.cell.ids = c("nht3", "nht5","nht9"), project = "nHT",merge.data = TRUE)

scRNAsub <- pbmc
scRNAsub$seurat_clusters[scRNAsub$seurat_clusters%in%c(0,1,2,3)]<-"0_1_2_3"
data <- as(as.matrix(scRNAsub@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNAsub@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
# Now, make a new CellDataSet using the RNA counts
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())   ###转为了RPS，分布模型改为negbinomial.size()
rm(data)
mycds <- estimateSizeFactors(mycds)
##使用nTE和mTE中HT和非HT的高比例亚群的前50特征基因
##ht的nTE：0,2；ht的mTE：3；nht的nTE：1；nht的mTE：0,1,2；
Idents(scRNAsub)<-"seurat_clusters"
diff.genes<- FindAllMarkers(object = scRNAsub,
                            only.pos = FALSE,
                            min.pct = 0.25,
                            logfc.threshold = 0)
cluster_diff<-diff.genes
cluster_diff <- subset(cluster_diff,p_val_adj<0.05&abs(avg_log2FC)>0.5)$gene
rm(scRNAsub)


mycds <- estimateDispersions(mycds, cores=8, relative_expr = TRUE)
mycds <- setOrderingFilter(mycds, cluster_diff)
#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree"')

#排序
mycds <- orderCells(mycds)
#State轨迹分布图
plot1 <- plot_cell_trajectory(mycds, color_by = "State")+
  scale_colour_manual(
    values = mypal[1:9]
    # aesthetics = c("colour", "fill")
  )
ggsave("State.pdf", plot = plot1, width = 6, height = 5)
##Cluster轨迹分布图

cols_1<-jjAnno::useMyCol("stallion",n=20)
plot2 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters")+
  scale_colour_manual(
    values = cols_1
    # aesthetics = c("colour", "fill")
  )
ggsave("seurat_clusters.pdf", plot = plot2, width = 6, height = 5)
##Pseudotime轨迹图
plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
ggsave("Pseudotime.pdf", plot = plot3, width = 6, height = 5)

p2 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters")+NoLegend()+
  scale_colour_manual(
    values = cols_1
    # aesthetics = c("colour", "fill")
  ) + facet_wrap(~seurat_clusters, nrow =2)

ggsave("trajectory_facet.pdf", plot = p2, width = 6, height =6)


####选择用差异基因
#cluster差异基因
#diff.genes <- read.csv('subcluster/diff_genes_wilcox.csv')
sig_diff.genes <- c.TopMarkers[which(c.TopMarkers$cluster%in%c(3)),"gene"]#用我们关注的亚群的marker基因
sig_diff.genes <- unique(as.character(sig_diff.genes$gene))
diff_test <- differentialGeneTest(mycds[sig_diff.genes,], cores = 1, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test, qval < 0.05))
p1 = plot_pseudotime_heatmap(mycds[sig_gene_names,], num_clusters=5,
                             show_rownames=T, return_heatmap=T)
ggsave("pseudotime_heatmap1.pdf", plot = p1, width = 5, height = 8)


mycds_sub <- mycds[sig_gene_names,]
beam_res <- BEAM(mycds_sub, branch_point = 1, cores = 8)
beam_res <- beam_res[order(beam_res$qval),]
beam_res <- beam_res[,c("gene_short_name", "pval", "qval")]
mycds_sub_beam <- mycds_sub[row.names(subset(beam_res, qval < 0.05)),]
pdf("pseudotime_point1.pdf",width = 5,height = 8)
plot_genes_branched_heatmap(mycds_sub_beam,  branch_point = 1, num_clusters =6,show_rownames = T)
dev.off()

mycds_sub <- mycds[sig_gene_names,]
beam_res <- BEAM(mycds_sub, branch_point = 2, cores = 8)
beam_res <- beam_res[order(beam_res$qval),]
beam_res <- beam_res[,c("gene_short_name", "pval", "qval")]
mycds_sub_beam <- mycds_sub[row.names(subset(beam_res, qval < 0.05)),]
pdf("pseudotime_point2.pdf",width = 5,height = 8)
plot_genes_branched_heatmap(mycds_sub_beam,  branch_point = 2, num_clusters = 2, show_rownames = T)
dev.off()

mycds_sub <- mycds[sig_gene_names,]
beam_res <- BEAM(mycds_sub, branch_point = 3, cores = 8)
beam_res <- beam_res[order(beam_res$qval),]
beam_res <- beam_res[,c("gene_short_name", "pval", "qval")]
mycds_sub_beam <- mycds_sub[row.names(subset(beam_res, qval < 0.05)),]
pdf("pseudotime_point3.pdf",width = 5,height = 8)
plot_genes_branched_heatmap(mycds_sub_beam,  branch_point = 3, num_clusters = 2, show_rownames = T)
dev.off()

mycds_sub <- mycds[sig_gene_names,]
beam_res <- BEAM(mycds_sub, branch_point = 4, cores = 8)
beam_res <- beam_res[order(beam_res$qval),]
beam_res <- beam_res[,c("gene_short_name", "pval", "qval")]
mycds_sub_beam <- mycds_sub[row.names(subset(beam_res, qval < 0.05)),]
pdf("pseudotime_point4.pdf",width = 5,height = 8)
plot_genes_branched_heatmap(mycds_sub_beam,  branch_point =4, num_clusters =2, show_rownames = T)
dev.off()

a <- subset(scRNAsub.T@meta.data, seurat_clusters==4)
a<- row.names(a)
write.table(a,"Fibroblasts_4.txt",col.names = F,row.names = F,sep = "\t",quote = F)
b<-c.TopMarkers[which(c.TopMarkers$cluster==4),"gene"]
b<-b$gene
write.table(b,"Fibroblasts_4_marker.txt",col.names = F,row.names = F,sep = "\t",quote = F)
