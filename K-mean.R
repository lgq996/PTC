library(infercnv)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")
library(ggsci)
mypal<-pal_npg("nrc")(10)

infercnv_obj = readRDS("run.final.infercnv_obj")
expr <- infercnv_obj@expr.data
normal_loc <- infercnv_obj@reference_grouped_cell_indices
normal_loc <- normal_loc$`para-tumors(non-malignant)`
test_loc <- infercnv_obj@observation_grouped_cell_indices
test_loc <- c(test_loc$`distant-metastases`,test_loc$`Lymph-Node`,test_loc$`primary-tumors`)

anno.df=data.frame(
  CB=c(colnames(expr)[normal_loc],colnames(expr)[test_loc]),
  class=c(rep("normal",length(normal_loc)),rep("test",length(test_loc)))
)
head(anno.df)

gn <- rownames(expr)
geneFile <- read.table("geneLocate.txt",header = F,sep = "\t",stringsAsFactors = F)
rownames(geneFile)=geneFile$V1
sub_geneFile <-  geneFile[intersect(gn,geneFile$V1),]
expr=expr[intersect(gn,geneFile$V1),]
head(sub_geneFile,4)
expr[1:4,1:4]

#聚类，10类，提取结果
set.seed(20210418)
for(i in 2){
  kmeans.result <- kmeans(t(expr), i)
  kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
  kmeans_df$CB=rownames(kmeans_df)
  kmeans_df=kmeans_df%>%inner_join(anno.df,by="CB") #合并
  kmeans_df_s=arrange(kmeans_df,kmeans_class) #排序
  rownames(kmeans_df_s)=kmeans_df_s$CB
  kmeans_df_s$CB=NULL
  kmeans_df_s$kmeans_class=as.factor(kmeans_df_s$kmeans_class) #将kmeans_class转换为因子，作为热图的一个注释，最终在热图上就会按照1:7排序
  head(kmeans_df_s)
  
  #定义热图的注释，及配色
  top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
  color_v=colorRampPalette(colors = c(mypal))(i) #类别数
  names(color_v)=as.character(1:i)
  left_anno <- rowAnnotation(df = kmeans_df_s,col=list(class=c("test"="red","normal" = "blue"),kmeans_class=color_v))
  
  #下面是绘图
  pdf(paste0("try_",i,".pdf"),width = 15,height = 10)
  ht = Heatmap(t(expr)[rownames(kmeans_df_s),], #绘图数据的CB顺序和注释CB顺序保持一致
               col = colorRamp2(c(0.4,1,1.6), c("#377EB8","#F0F0F0","#E41A1C")), #如果是10x的数据，这里的刻度会有所变化
               cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
               column_split = factor(sub_geneFile$V2, paste("chr",1:22,sep = "")), #这一步可以控制染色体顺序，即使你的基因排序文件顺序是错的
               column_gap = unit(2, "mm"),
               
               heatmap_legend_param = list(title = "Modified expression",direction = "vertical",title_position = "leftcenter-rot",at=c(0.4,1,1.6),legend_height = unit(3, "cm")),
               
               top_annotation = top_anno,left_annotation = left_anno, #添加注释
               row_title = NULL,column_title = NULL)
  draw(ht, heatmap_legend_side = "right")
  dev.off()
  
  #每一类对应的CB保存在kmeans_df_s数据框中
  write.table(kmeans_df_s, file =paste0("kmeans_df_s",i,".txt"), quote = FALSE, sep = '\t', row.names = T, col.names = T)
  
}
