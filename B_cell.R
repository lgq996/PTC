load(file = "B.Rdata")
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

memory.limit(size = 100000)

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
scRNAsub.T <- FindClusters(object = scRNAsub.T, resolution = 0.8)###粗聚类0.3，细聚类0.8                 
scRNAsub.T <- RunTSNE(object = scRNAsub.T, dims = 1:pcSelect)                      
pdf(file="06.TSNE_Bcell.pdf",width=8,height=6)
DimPlot(object = scRNAsub.T, label = TRUE,pt.size = 1,label.size = 10,label.box = T,repel =T,group.by = "seurat_clusters",
        cols = colorRampPalette(colors = c(mypal))(20))+NoLegend()+
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
features<-c("MS4A1","IGHD","CD27","NEIL1","MZB1","CD38")

p1<-DotPlot(scRNAsub.T,features = features,cols = colorRampPalette(colors = c(mypal))(20),group.by = "cellType_B")
pdot<-p1$data
pdot$features.plot<-factor(pdot$features.plot,levels = unique(pdot$features.plot))

pdot<-ggplot(pdot,aes(x=features.plot,y=id))+
  geom_point(aes(fill=avg.exp,size=pct.exp),color='black',shape=21)+
  theme_bw(base_size = 12)+
  xlab("")+ylab("")+
  scale_fill_gradient2(low = '#4DBBD5FF',mid = 'white',high = '#E64B35FF',midpoint = 1,name='Mean expression')+
  scale_size(range = c(0.5,10),name='Percentage')+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        aspect.ratio = 0.5,
        plot.margin = margin(t=1,r=1,b=1,l=1,unit = 'cm'),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,face = 'bold'),
        axis.text.y = element_text(angle = 0,hjust = 1,vjust = 0.5,face = 'bold'))+
  coord_cartesian(clip = 'off')
pdot_test<-pdot+
  scale_x_discrete(position = 'top')+
  theme(axis.text.x = element_text(angle = 90,hjust = 0))
ggsave("cell_bubble_B_celltype.pdf",width = 10,height =8)
###绘制亚群间差异基因,主要目的是进一步确定亚群的细胞类型
scRNAsub.T@meta.data$cellType_B<-ifelse(scRNAsub.T@meta.data$seurat_clusters%in%c(0,1,2,7,9,10,12,14),"Native_B",
                                        ifelse(scRNAsub.T@meta.data$seurat_clusters%in%c(3,4,6,8,11,13,17),"Memory_B",
                                               ifelse(scRNAsub.T@meta.data$seurat_clusters%in%c(5,15,19),"Intermediate_B",
                                                      ifelse(scRNAsub.T@meta.data$seurat_clusters%in%c(16),"Plasma_B","Germinal_center_B"))))
#先画一个，获取图例
pdf(file="06.TSNE_Bcell_cellType.pdf",width=8,height=6)
DimPlot(object = scRNAsub.T, label = T,pt.size = 1,label.size =6,label.box = T,repel =T,group.by = "cellType_B",
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
#按照细胞类型做差异分析
Idents(scRNAsub.T)="cellType_B"
skcm.markers_celltype.B <- FindAllMarkers(object = scRNAsub.T,
                                 only.pos = FALSE,
                                 min.pct = 0.25,
                                 logfc.threshold = 0)
top = 50 # can be adjusted as needed
TopMarkers <- skcm.markers_celltype.B %>% filter(p_val_adj < 0.01 & abs(avg_log2FC)>0.5) %>% group_by(cluster) %>%
  top_n(n = top, wt = avg_log2FC)
###画细胞比例堆积条形图
#devtools::install_github("caleblareau/BuenColors")
library(reshape2)
library(BuenColors)
jdb_palette("lawhoops")
bar<-scRNAsub.T@meta.data
bar<-bar[,c(8,9,13)]
bar<-reshape2::dcast(bar,HT+group~cellType_B,fill = 0)
bar<-melt(bar)
#bar$Group<-factor(bar$Group,levels = c("Pre","Post"))
ggplot(data = bar,aes(x=HT,y=value))+
  geom_bar(aes(fill=variable),stat = 'identity',position = position_fill(),width = 0.5)+
  theme_bw()+scale_fill_manual(values =c("Native_B"=jdb_palette("lawhoops")[1],"Memory_B"=jdb_palette("lawhoops")[3],
                                         "Plasma_B"=jdb_palette("lawhoops")[5],"Intermediate_B"=jdb_palette("lawhoops")[8],
                                         "Germinal_center_B"=jdb_palette("lawhoops")[7]))+
  #facet_grid(~group,scales = "free")+
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
ggsave("bar_B.pdf",width =6,height = 6)
library(rstatix)
library(SeuratDisk)
cellinfo <- Seurat::FetchData(scRNAsub.T, vars = c("HT", "group","cellType_B"))
#cellinfo<-cellinfo[which(cellinfo$group=="Normal"),]
cell.col <- setNames(object = c(mypal[3], mypal[9]),
                     nm = c("HT", "Non-HT"))

# 构建seurat_clusters和Group的列联表
tbl <- table(cellinfo$cellType_B, cellinfo$HT)
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
ggsave(filename = "Bbar_sig.pdf", width = 8, height = 6)


###细胞通讯分析
library(Seurat)
library(tidyverse)
library(iTALK)

###HT样本的T(患者1,2,8)
ht1<-subset(x=scRNAsub.T,orig.ident%in%c("PTC1P","PTC1T"))#1和2的细胞数太少，此处用8的
ht2<-subset(x=scRNAsub.T,orig.ident%in%c("PTC2LeftLN","PTC2P","PTC2T"))
ht8<-subset(x=scRNAsub.T,orig.ident%in%c("PTC8P","PTC8T"))
###非HT样本的T(患者3,5,9)
nht3<-subset(x=scRNAsub.T,orig.ident%in%c("PTC3LeftLN","PTC3P","PTC3RightLN","PTC3T"))
nht5<-subset(x=scRNAsub.T,orig.ident%in%c("PTC5P","PTC5RightLN","PTC5T"))
nht9<-subset(x=scRNAsub.T,orig.ident%in%c("PTC9P","PTC9T"))

##==细胞通讯关系概览==##
cell.meta <- subset(ht8@meta.data, select=c("orig.ident","cellType_B"))
names(cell.meta) <- c("compare_group", "cell_type")
cell.expr <- data.frame(t(as.matrix(ht8@assays$RNA@counts)), check.names=F)
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

paitent<-"HT_PTC8"
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
  setwd("E:/lgq/重新分析/B_cell")
}
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
Idents(scRNAsub)<-"cellType_B"
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
plot2 <- plot_cell_trajectory(mycds, color_by = "cellType_B")+
  scale_colour_manual(
    values = cols_1
    # aesthetics = c("colour", "fill")
  )
ggsave("celltype_B.pdf", plot = plot2, width = 6, height = 5)
##Pseudotime轨迹图
plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
ggsave("Pseudotime.pdf", plot = plot3, width = 6, height = 5)

p2 <- plot_cell_trajectory(mycds, color_by = "cellType_B")+NoLegend()+
  scale_colour_manual(
    values = cols_1
    # aesthetics = c("colour", "fill")
  ) + facet_wrap(~cellType_B, nrow =2)

ggsave("trajectory_facet.pdf", plot = p2, width = 8, height =6)


####选择用差异基因
#cluster差异基因
#diff.genes <- read.csv('subcluster/diff_genes_wilcox.csv')
sig_diff.genes <- TopMarkers[which(TopMarkers$cluster%in%c("Germinal_center_B")),"gene"]#用我们关注的亚群的marker基因
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
plot_genes_branched_heatmap(mycds_sub_beam,  branch_point = 2, num_clusters = 3, show_rownames = T)
dev.off()

mycds_sub <- mycds[sig_gene_names,]
beam_res <- BEAM(mycds_sub, branch_point = 3, cores = 8)
beam_res <- beam_res[order(beam_res$qval),]
beam_res <- beam_res[,c("gene_short_name", "pval", "qval")]
mycds_sub_beam <- mycds_sub[row.names(subset(beam_res, qval < 0.05)),]
pdf("pseudotime_point3.pdf",width = 5,height = 8)
plot_genes_branched_heatmap(mycds_sub_beam,  branch_point = 3, num_clusters = 3, show_rownames = T)
dev.off()

mycds_sub <- mycds[sig_gene_names,]
beam_res <- BEAM(mycds_sub, branch_point = 4, cores = 8)
beam_res <- beam_res[order(beam_res$qval),]
beam_res <- beam_res[,c("gene_short_name", "pval", "qval")]
mycds_sub_beam <- mycds_sub[row.names(subset(beam_res, qval < 0.05)),]
pdf("pseudotime_point4.pdf",width = 5,height = 8)
plot_genes_branched_heatmap(mycds_sub_beam,  branch_point =4, num_clusters =3, show_rownames = T)
dev.off()

a <- subset(scRNAsub.T@meta.data, cellType_B=="Intermediate_B")
a<- row.names(a)
write.table(a,"Intermediate_B.txt",col.names = F,row.names = F,sep = "\t",quote = F)
b<-TopMarkers[which(TopMarkers$cluster=="Intermediate_B"),"gene"]
b<-b$gene
write.table(b,"Intermediate_B_marker.txt",col.names = F,row.names = F,sep = "\t",quote = F)
