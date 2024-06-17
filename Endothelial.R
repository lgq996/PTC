load(file = "Endothelial.Rdata")
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
scRNAsub.T <- FindClusters(object = scRNAsub.T, resolution = 0.7)###粗聚类0.3，细聚类0.8                 
scRNAsub.T <- RunTSNE(object = scRNAsub.T, dims = 1:pcSelect)                      
pdf(file="06.TSNE_Ecell.pdf",width=8,height=6)
DimPlot(object = scRNAsub.T, label = TRUE,pt.size = 1,label.size = 10,label.box = T,repel =T,group.by = "seurat_clusters",
        cols = colorRampPalette(colors = c(mypal))(15))+NoLegend()+
  labs(x="tSNE1",y="tSNE2",title = "")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 15,face = "bold"), 
        axis.text.y = element_text(angle = 0, size = 15,face = "bold"),
        axis.title.y = element_text(angle = 90, size = 20,face = "bold"),
        axis.title.x = element_text(angle = 0, size = 20,face = "bold"),
        legend.title = element_text(size = 20,face =  "bold"),
        legend.text = element_text(size = 20,face = "bold"),
        panel.border = element_rect(fill = NA,colour = "black",size = 2,linetype = "solid"))
dev.off()
###通过富集分析进一步确定亚型
Idents(scRNAsub.T)<-"seurat_clusters"
markers.CAF.cluster<- FindAllMarkers(object = scRNAsub.T,
                                  only.pos = FALSE,
                                  min.pct = 0.25,
                                  logfc.threshold = 0)
top = 50 # can be adjusted as needed
TopMarkers <- markers.CAF.cluster %>% filter(p_val_adj < 0.01 & abs(avg_log2FC)>0.5) %>% group_by(cluster) %>%
  top_n(n = top, wt = avg_log2FC)
###根据每个cluster的marker基因可以做富集分析，来确定亚群的生物学作用
###绘制热图
library(ClusterGVis)
library(org.Hs.eg.db)
# prepare data from seurat object
st.data <- prepareDataFromscRNA(object = scRNAsub.T,
                                diffData = TopMarkers,
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
col<-c(rep(jjAnno::useMyCol("stallion",n = 15),each = 5))
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
           cluster.order = c(1:15),
           go.col = col,
           add.bar = T)
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
  theme_bw()+scale_fill_manual(values =c(colorRampPalette(colors = c(jdb_palette("lawhoops")))(15)))+
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
ggsave("bar_E_all.pdf",width =6,height = 6)
library(rstatix)
library(SeuratDisk)
cellinfo <- Seurat::FetchData(scRNAsub.T, vars = c("HT", "group","seurat_clusters"))
#cellinfo<-cellinfo[which(cellinfo$group=="metastasis"),]
cell.col <- setNames(object = c(mypal[9], mypal[10]),
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
ggsave(filename = "Ebar_sig.pdf", width = 8, height = 6)
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
scRNAsub$seurat_clusters[scRNAsub$seurat_clusters%in%c(1,3,4,5)]<-"1_3_4_5"
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
  ) + facet_wrap(~seurat_clusters, nrow =3)

ggsave("trajectory_facet.pdf", plot = p2, width = 6, height =8)


####选择用差异基因
#cluster差异基因
#diff.genes <- read.csv('subcluster/diff_genes_wilcox.csv')
sig_diff.genes <- TopMarkers[which(TopMarkers$cluster%in%c(5)),"gene"]#用我们关注的亚群的marker基因
sig_diff.genes <- unique(as.character(sig_diff.genes$gene))
diff_test <- differentialGeneTest(mycds[sig_diff.genes,], cores = 1, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test, qval < 0.05))
p1 = plot_pseudotime_heatmap(mycds[sig_gene_names,], num_clusters=4,
                             show_rownames=T, return_heatmap=T)
ggsave("pseudotime_heatmap1.pdf", plot = p1, width = 5, height = 8)


mycds_sub <- mycds[sig_gene_names,]
beam_res <- BEAM(mycds_sub, branch_point = 1, cores = 8)
beam_res <- beam_res[order(beam_res$qval),]
beam_res <- beam_res[,c("gene_short_name", "pval", "qval")]
mycds_sub_beam <- mycds_sub[row.names(subset(beam_res, qval < 0.05)),]
pdf("pseudotime_point1.pdf",width = 5,height = 8)
plot_genes_branched_heatmap(mycds_sub_beam,  branch_point = 1, num_clusters =3,show_rownames = T)
dev.off()

mycds_sub <- mycds[sig_gene_names,]
beam_res <- BEAM(mycds_sub, branch_point = 2, cores = 8)
beam_res <- beam_res[order(beam_res$qval),]
beam_res <- beam_res[,c("gene_short_name", "pval", "qval")]
mycds_sub_beam <- mycds_sub[row.names(subset(beam_res, qval < 0.05)),]
pdf("pseudotime_point2.pdf",width = 5,height = 8)
plot_genes_branched_heatmap(mycds_sub_beam,  branch_point = 2, num_clusters = 3, show_rownames = T)
dev.off()

mycds_sub <- mycds[sig_gene_names,]
beam_res <- BEAM(mycds_sub, branch_point = 3, cores = 8)
beam_res <- beam_res[order(beam_res$qval),]
beam_res <- beam_res[,c("gene_short_name", "pval", "qval")]
mycds_sub_beam <- mycds_sub[row.names(subset(beam_res, qval < 0.05)),]
pdf("pseudotime_point3.pdf",width = 5,height = 8)
plot_genes_branched_heatmap(mycds_sub_beam,  branch_point = 3, num_clusters = 4, show_rownames = T)
dev.off()

mycds_sub <- mycds[sig_gene_names,]
beam_res <- BEAM(mycds_sub, branch_point = 4, cores = 8)
beam_res <- beam_res[order(beam_res$qval),]
beam_res <- beam_res[,c("gene_short_name", "pval", "qval")]
mycds_sub_beam <- mycds_sub[row.names(subset(beam_res, qval < 0.05)),]
pdf("pseudotime_point4.pdf",width = 5,height = 8)
plot_genes_branched_heatmap(mycds_sub_beam,  branch_point =4, num_clusters =2, show_rownames = T)
dev.off()

a <- subset(scRNAsub.T@meta.data, seurat_clusters==13)
a<- row.names(a)
write.table(a,"Endothelial_13.txt",col.names = F,row.names = F,sep = "\t",quote = F)
b<-TopMarkers[which(TopMarkers$cluster==13),"gene"]
b<-b$gene
write.table(b,"Endothelial_13_marker.txt",col.names = F,row.names = F,sep = "\t",quote = F)
