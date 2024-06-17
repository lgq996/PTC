load("TE_inferCNV.Rdata")
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
scRNAsub.T$inferCNV[scRNAsub.T$inferCNV=="para-tumors-tumors"]<-"primary-tumors"
###提取恶性甲状腺细胞
subT<-subset(x=scRNAsub.T,inferCNV=="para-tumors",invert = TRUE)
subT <- FindVariableFeatures(object = subT, selection.method = "vst", nfeatures = 5000)

subT=ScaleData(subT)                    
subT=RunPCA(object= subT,npcs = 20,pc.genes=VariableFeatures(object = subT))    

subT <- JackStraw(object = subT, num.replicate = 100)
subT <- ScoreJackStraw(object = subT, dims = 1:20)
JackStrawPlot(object = subT, dims = 1:20,
              cols = colorRampPalette(colors = c(mypal))(20))
###################################06.TSNE??????????marker????###################################
##TSNE????????
pcSelect=20
subT <- FindNeighbors(object = subT, dims = 1:pcSelect)
#先粗略聚类
subT <- FindClusters(object = subT, resolution = 0.5)###粗聚类0.3，细聚类0.8                 
subT <- RunTSNE(object = subT, dims = 1:pcSelect)                      
pdf(file="06.TSNE_TEcell.pdf",width=8,height=6)
DimPlot(object = subT, label = TRUE,pt.size = 1,label.size = 10,label.box = T,repel =T,group.by = "seurat_clusters",
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
###绘制亚群间差异基因,主要目的是进一步确定亚群的细胞类型
skcm.markers.TE <- FindAllMarkers(object = subT,
                                 only.pos = FALSE,
                                 min.pct = 0.25,
                                 logfc.threshold = 0)
top = 50 # can be adjusted as needed
TopMarkers <- skcm.markers.TE %>% filter(p_val_adj < 0.01 & abs(avg_log2FC)>0.5) %>% group_by(cluster) %>%
  top_n(n = top, wt = avg_log2FC)
###根据每个cluster的marker基因可以做富集分析，来确定亚群的生物学作用
###绘制热图
library(ClusterGVis)
library(org.Hs.eg.db)
# prepare data from seurat object
st.data <- prepareDataFromscRNA(object = subT,
                                diffData = TopMarkers,
                                showAverage = F, keep.uniqGene = FALSE,
                                sep = "_")
library(R.utils)

R.utils::setOption("clusterProfiler.download.method",'auto')
# enrich for clusters
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "MF",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 5,
                        seed = 5201314)
col<-c(rep(jjAnno::useMyCol("stallion",n = 14),each = 5))
col<-col[-c(64,65)]
# add GO annotation
pdf('GO_MF.pdf',height = 16,width = 20,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           column_title_rot = 90,
           markGenes = NULL,#基因数太多，此处不显示
           markGenes.side = "left",
           annoTerm.data = enrich,
           show_column_names = F,
           line.side = "left",
           show_row_dend = F,
           cluster.order = c(1:14),
           go.col = col,
           add.bar = T)
dev.off()
###确定在HT和非HT之间各个亚群的占比
#按照细胞类型做差异分析
###画细胞比例堆积条形图
#devtools::install_github("caleblareau/BuenColors")
library(reshape2)
library(BuenColors)
jdb_palette("lawhoops")
bar<-subT@meta.data
#癌旁组织的肿瘤就是原发肿瘤
#bar$inferCNV<-ifelse(bar$inferCNV=="para-tumors","para-tumors","tumors")

bar<-bar[,c(8,10,6)]
bar<-reshape2::dcast(bar,HT~seurat_clusters,fill = 0)
bar<-melt(bar)
#bar$Group<-factor(bar$Group,levels = c("Pre","Post"))
ggplot(data = bar,aes(x=HT,y=value))+
  geom_bar(aes(fill=variable),stat = 'identity',position = position_fill(),width = 0.5)+
  theme_bw()+scale_fill_manual(values =c(colorRampPalette(colors = c(jdb_palette("lawhoops")))(15)))+
  #facet_grid(~inferCNV,scales = "free")+
  labs(x="",y="",title = "")+labs(fill="mTE_clusters")+
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
ggsave("bar_1.pdf",width =6,height = 6)
library(rstatix)
library(SeuratDisk)
cellinfo <- Seurat::FetchData(subT, vars = c("HT", "seurat_clusters"))
#cellinfo<-cellinfo[which(cellinfo$group=="Normal"),]
cell.col <- setNames(object = c(mypal[8], mypal[9]),
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
plot.data$Var1<-factor(plot.data$Var1,levels = c("mTE0","mTE1","mTE2","mTE3","mTE4","mTE5","mTE6",
                                                 "mTE7","mTE8","mTE9","mTE10","mTE11","mTE12","mTE13"))
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
        axis.text = element_text(size = 15,face = "bold"),              # 调整坐标轴刻度字号大小
        axis.text.x = element_text(face = "bold",angle = 45,hjust = 1),
        axis.title = element_text(size = 20,face = "bold"),
        legend.title = element_text(size = 20,face =  "bold"),          # 调整坐标轴标题字号大小
        legend.text = element_text(size = 10,face = "bold"),
  )                           
p2
ggsave(filename = "mbar_sig.pdf", width = 10, height = 6)
save(subT,file = "mTE.Rdata")

###提取恶性甲状腺细胞
subN<-subset(x=scRNAsub.T,inferCNV=="para-tumors")
subN <- FindVariableFeatures(object = subN, selection.method = "vst", nfeatures = 5000)

subN=ScaleData(subN)                    
subN=RunPCA(object= subN,npcs = 20,pc.genes=VariableFeatures(object = subN))    

subN <- JackStraw(object = subN, num.replicate = 100)
subN <- ScoreJackStraw(object = subN, dims = 1:20)
JackStrawPlot(object = subN, dims = 1:20,
              cols = colorRampPalette(colors = c(mypal))(20))
###################################06.TSNE??????????marker????###################################
##TSNE????????
pcSelect=20
subN <- FindNeighbors(object = subN, dims = 1:pcSelect)
#先粗略聚类
subN <- FindClusters(object = subN, resolution = 0.5)###粗聚类0.3，细聚类0.8                 
subN <- RunTSNE(object = subN, dims = 1:pcSelect)                      
pdf(file="06.TSNE_nTEcell.pdf",width=8,height=6)
DimPlot(object = subN, label = TRUE,pt.size = 1,label.size = 10,label.box = T,repel =T,group.by = "seurat_clusters",
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
###绘制亚群间差异基因,主要目的是进一步确定亚群的细胞类型
skcm.markers.nTE <- FindAllMarkers(object = subN,
                                  only.pos = FALSE,
                                  min.pct = 0.25,
                                  logfc.threshold = 0)
top = 50 # can be adjusted as needed
nTopMarkers <- skcm.markers.nTE %>% filter(p_val_adj < 0.01 & abs(avg_log2FC)>0.5) %>% group_by(cluster) %>%
  top_n(n = top, wt = avg_log2FC)
###根据每个cluster的marker基因可以做富集分析，来确定亚群的生物学作用
###绘制热图
library(ClusterGVis)
library(org.Hs.eg.db)
# prepare data from seurat object
nst.data <- prepareDataFromscRNA(object = subN,
                                diffData = nTopMarkers,
                                showAverage = F, keep.uniqGene = FALSE,
                                sep = "_")
library(R.utils)

R.utils::setOption("clusterProfiler.download.method",'auto')
# enrich for clusters
enrich <- enrichCluster(object = nst.data,
                        OrgDb = org.Hs.eg.db,
                        type = "MF",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 5,
                        seed = 5201314)
col<-c(rep(jjAnno::useMyCol("stallion",n = 11),each = 5))
col<-col[-c(64,65)]
# add GO annotation
pdf('nGO_MF.pdf',height = 16,width = 18,onefile = F)
visCluster(object = nst.data,
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
###确定在HT和非HT之间各个亚群的占比
#按照细胞类型做差异分析
###画细胞比例堆积条形图
#devtools::install_github("caleblareau/BuenColors")
library(reshape2)
library(BuenColors)
jdb_palette("lawhoops")
bar<-subN@meta.data
#癌旁组织的肿瘤就是原发肿瘤
#bar$inferCNV<-ifelse(bar$inferCNV=="para-tumors","para-tumors","tumors")

bar<-bar[,c(8,6)]
bar<-reshape2::dcast(bar,HT~seurat_clusters,fill = 0)
bar<-melt(bar)
#bar$Group<-factor(bar$Group,levels = c("Pre","Post"))
ggplot(data = bar,aes(x=HT,y=value))+
  geom_bar(aes(fill=variable),stat = 'identity',position = position_fill(),width = 0.5)+
  theme_bw()+scale_fill_manual(values =c(colorRampPalette(colors = c(jdb_palette("lawhoops")))(15)))+
  #facet_grid(~inferCNV,scales = "free")+
  labs(x="",y="",title = "")+labs(fill="nTE_clusters")+
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
ggsave("nbar_1.pdf",width =6,height = 6)
library(rstatix)
library(SeuratDisk)
cellinfo <- Seurat::FetchData(subN, vars = c("HT", "seurat_clusters"))
#cellinfo<-cellinfo[which(cellinfo$group=="Normal"),]
cell.col <- setNames(object = c(mypal[3], mypal[4]),
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
plot.data$Var1<-factor(plot.data$Var1,levels = c("nTE0","nTE1","nTE2","nTE3","nTE4","nTE5","nTE6",
                                                 "nTE7","nTE8","nTE9","nTE10"))
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
        axis.text = element_text(size = 15,face = "bold"),              # 调整坐标轴刻度字号大小
        axis.text.x = element_text(face = "bold",angle = 45,hjust = 1),
        axis.title = element_text(size = 20,face = "bold"),
        legend.title = element_text(size = 20,face =  "bold"),          # 调整坐标轴标题字号大小
        legend.text = element_text(size = 10,face = "bold"),
  )                           
p2
ggsave(filename = "nbar_sig.pdf", width = 10, height = 6)
save(subN,file = "nTE.Rdata")
###分别分析HT和非HT样本的nTE细胞和mTE细胞之间的通讯以及轨迹分化
###在nTE和mTE之间是如何进行细胞通讯级轨迹分化的，此处将不会考虑远端，因为差异的亚群主要分布在原发肿瘤
###先合并N和T
#分别标记两类亚群
subN$seurat_clusters<-paste0("nTE",subN$seurat_clusters)
subT$seurat_clusters<-paste0("mTE",subT$seurat_clusters)

pbmc <- merge(subN, y = subT, add.cell.ids = c("subN", "subT"), project = "ALL",merge.data = TRUE)
###首先分析HT样本
###HT样本的nTE(患者1,2,8)
ht<-subset(x=pbmc,orig.ident%in%c("PTC8P","PTC8T"))#1和2的细胞数太少，此处用8的
nht5<-subset(x=pbmc,orig.ident%in%c("PTC5P","PTC5RightLN","PTC5T"))
nht9<-subset(x=pbmc,orig.ident%in%c("PTC9P","PTC9T"))

table(ht$seurat_clusters)
#亚群mTE3和mTE10为有利于甲状腺癌的亚群，重点关注哪些亚群与这两个存在通讯及哪些亚群会分化为这两个亚群
#通讯分析
library(Seurat)
library(tidyverse)
library(iTALK)
library(dplyr)
##准备iTALK输入文件
# 输入文件格式为数据框，行为细胞名称，列为基因名和细胞注释信息，基因名命名的列值为表达数据。
# 细胞注释信息必选“cell_type”，值为细胞类型注释信息；可选“compare_group”，值为细胞的样本分组信息。

##==细胞通讯关系概览==##
cell.meta <- subset(nht9@meta.data, select=c("orig.ident","seurat_clusters"))
names(cell.meta) <- c("compare_group", "cell_type")
cell.expr <- data.frame(t(as.matrix(nht9@assays$RNA@counts)), check.names=F)
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

paitent<-"nHT_PTC9"
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
  setwd("E:/lgq/重新分析/TE/亚群聚类")
}
###轨迹分析
library(Seurat)
library(dplyr)
library(monocle)
library(ggsci)
#mypal<-pal_jco("default", alpha = 0.6)(9)
scRNAsub <- nht9

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
    values = mypal[1:7]
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
  ) + facet_wrap(~seurat_clusters, nrow =4)

ggsave("trajectory_facet.pdf", plot = p2, width = 8, height =6)


####选择用差异基因
#cluster差异基因
#diff.genes <- read.csv('subcluster/diff_genes_wilcox.csv')
sig_diff.genes <- nTopMarkers[which(nTopMarkers$cluster%in%c(1)),"gene"]#用我们关注的亚群的marker基因
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
plot_genes_branched_heatmap(mycds_sub_beam,  branch_point = 1, num_clusters =2,show_rownames = T)
dev.off()

mycds_sub <- mycds[sig_gene_names,]
beam_res <- BEAM(mycds_sub, branch_point = 2, cores = 8)
beam_res <- beam_res[order(beam_res$qval),]
beam_res <- beam_res[,c("gene_short_name", "pval", "qval")]
mycds_sub_beam <- mycds_sub[row.names(subset(beam_res, qval < 0.05)),]
pdf("pseudotime_point2.pdf",width = 5,height = 8)
plot_genes_branched_heatmap(mycds_sub_beam,  branch_point = 2, num_clusters = 4, show_rownames = T)
dev.off()

mycds_sub <- mycds[sig_gene_names,]
beam_res <- BEAM(mycds_sub, branch_point = 3, cores = 8)
beam_res <- beam_res[order(beam_res$qval),]
beam_res <- beam_res[,c("gene_short_name", "pval", "qval")]
mycds_sub_beam <- mycds_sub[row.names(subset(beam_res, qval < 0.05)),]
pdf("pseudotime_point3.pdf",width = 5,height = 8)
plot_genes_branched_heatmap(mycds_sub_beam,  branch_point = 3, num_clusters = 5, show_rownames = T)
dev.off()


a <- subset(subN@meta.data, seurat_clusters=="nTE0")
a<- row.names(a)
write.table(a,"nTE0.txt",col.names = F,row.names = F,sep = "\t",quote = F)
b<-nTopMarkers[which(nTopMarkers$cluster==0),"gene"]
b<-b$gene
write.table(b,"nTE0_marker.txt",col.names = F,row.names = F,sep = "\t",quote = F)

a <- subset(subN@meta.data, seurat_clusters=="nTE2")
a<- row.names(a)
write.table(a,"nTE2.txt",col.names = F,row.names = F,sep = "\t",quote = F)
b<-nTopMarkers[which(nTopMarkers$cluster==2),"gene"]
b<-b$gene
write.table(b,"nTE2_marker.txt",col.names = F,row.names = F,sep = "\t",quote = F)

a <- subset(subT@meta.data, seurat_clusters=="mTE3")
a<- row.names(a)
write.table(a,"mTE3.txt",col.names = F,row.names = F,sep = "\t",quote = F)
b<-TopMarkers[which(TopMarkers$cluster==3),"gene"]
b<-b$gene
write.table(b,"mTE3_marker.txt",col.names = F,row.names = F,sep = "\t",quote = F)
