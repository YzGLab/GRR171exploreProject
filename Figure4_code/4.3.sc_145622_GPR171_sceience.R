###转录因子
subsubsce<-subsce[,subsce$assign.ident=="CD8T"]
###
write.csv(t(as.matrix(subsubsce@assays$RNA@counts)),file = "sce_exp.csv")
###
###
##
library(SCENIC)
# https://resources.aertslab.org/cistarget/
db='/home/data/gaoyuzhen/Projects/ImmuneProGPR/Data/cisTarget_database'
list.files(db)
#data(list="motifAnnotations_hgnc", package="RcisTarget")
# 保证cisTarget_databases 文件夹下面有下载好 的文件
motifAnnotations_hgnc <- motifAnnotations
scenicOptions <- initializeScenic(org="hgnc", 
                                  dbDir=db , nCores=16) 
#saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
#saveRDS(cellInfo, file="int/cellInfo.Rds")
data.input<-subsubsce@assays$RNA@data
#data.input<-data.input[,colnames(data.input) %in% fib.metadatas$Cell]
##
exprMat<-as.matrix(data.input) 
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered<-t(exprMat_filtered)

write.csv(exprMat_filtered,file = "sce_exp.csv")


rm(list = ls())  
library(SCENIC)
packageVersion("SCENIC")  
library(SCopeLoomR)
scenicLoomPath='/home/data/gaoyuzhen/Projects/ImmuneProGPR/Data/cisTarget_database/sample_SCENIC.loom'
loom <- open_loom(scenicLoomPath)
# Read information from loom file: 
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")##
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name="RegulonsAUC") 
###########################
Idents(subsubsce)<-subsubsce$group







#run Tcell SMC regulation
#sce<-readRDS("fibroblast_sce_annotation.RDS")
library(pheatmap) 
n=t(scale(t(getAUC(regulonAUC[,] )))) # 'scale'可以对log-ratio数值进行归一化
n[n>2]=2 
n[n< -2]= -2
n[1:4,1:4]
dim(n) 
ac=data.frame(group= as.character(Idents(subsubsce)))
n[1:4,1:4]
n=n[,colnames(n) %in% colnames(subsubsce)]
rownames(ac)=colnames(n) ##
lengths(regulons)
library(ggpubr)
#aucs<-getAUC(regulonAUC)
aucs<-getAUC(regulonAUC)
aucs<-t(aucs)
regulonss<-regulons
for(i in c(1:72)){
  names(regulonss)[i]<-gsub("\\+",paste0(lengths(regulons)[i],"g"),names(regulons)[i])
}
names(regulonss)
colnames(aucs)<-names(regulonss)
#aucs<-t(n)
#aucs<-data.frame(aucs)
aucs<-cbind(aucs,group= as.character(Idents(subsubsce)))
dim(aucs)
#write.csv(aucs,file="Bcell_TF_idente.csv")
#aucs<-read.csv("Bcell_TF_idente.csv",check.names = F)
head(aucs)

colnames(subsubsce)
aucs<-data.frame(aucs)
aucs[]<-apply(aucs,2,as.numeric)
subsubsce<-AddMetaData(subsubsce,aucs)

colnames(subsubsce@meta.data)[9:80]<-colnames(aucs)



subsubsce<-subsce[,subsce$assign.ident=="CD8T"]

subsubsce$group<-ifelse(subsubsce@assays$RNA@counts['GPR171',]>0,1,0)
subsubsce$group<-paste0(subsubsce$assign.ident,"_",subsubsce$group)


subsubsce<- RunUMAP(subsubsce, dims = 1:20)
CellDimPlot(
  srt = subsubsce, group.by =  "group",split.by = "group",
  cells.highlight = T,
  palcolor = c("black","red"),
  reduction = "UMAP", theme_use = "theme_blank"
)
#ETS2(44g)
FeatureDimPlot(
  srt =subsubsce, features = c("ETS2.44g."), split.by = "group",
  #cells.highlight = T,
  add_density = T,
  reduction = "UMAP", theme_use = "theme_blank"
)
data.frame(TF_heatmap)


library(ggplot2)
ht3 <- GroupHeatmap(
  srt =subsubsce,
  features =colnames(subsubsce@meta.data)[colnames(subsubsce@meta.data) %in% rownames(TF_heatmap)][-c(3,20,21,15:17)][-11],
  #palcolor =c("black","#E6AB02", "red"),
  #group_palette =NULL,
  #group_palcolor =c("#E6AB02", "red"),
  group.by =  "group",
  #heatmap_palette = "YlOrRd",
  heatmap_palette = "viridis",
  #cell_annotation = c( "cohort"),
  #cell_annotation_palette = c( "Paired"),
  show_row_names = TRUE, 
  row_names_side = "left",
  #add_violin = TRUE,
  #add_bg = TRUE,
  nlabel = 0,
  add_dot = TRUE, 
  #add_reticle = TRUE,
  border=TRUE
)
print(ht3$plot)
ggsave("Figure3/TF_T_GPR171.pdf",width = 4.9,height = 4.5)





#aucs<-data.frame(aucs)
#group_by(aucs, group) %>% summarize_each(funs(mean), colnames(aucs)[1:81])
aucss<-apply(aucs,2,as.numeric)
rownames(aucss)<-rownames(aucs)
aucsmenas<-aggregate(aucss[,c(1:72)],list(aucs[,73]),mean)#########根据最后确定的TF数量进行构建 heatmap
aucsmenas<-t(aucsmenas)
colnames(aucsmenas)<-aucsmenas[1,]
aucsmenas<-aucsmenas[-c(1),]
aucsmenas[1:5,1:4]
######
data<-aucs
outTab<-c()
for (i in colnames(aucs)[c(1:72)]){
  geneName=i
  data[,i]<-as.numeric(data[,i])
  rt=rbind(expression=data[,i],group=data[,"group"])
  rt=as.matrix(t(rt))
  rt<-data.frame(rt)
  #Means=rowMeans(data[,i])
  #if(is.numeric(Means) & Means>0){
  kTest<-kruskal.test(expression ~ group, data=rt)
  pvalue=kTest$p.value
  outTab=rbind(outTab,cbind(gene=i,pValue=pvalue))
  message(i,"is done now!!!")
  #}
}#
outTab<-data.frame(outTab)
#gene.cox.i[,11] <- p.adjust(gene.cox.i[,9], method = "fdr")
outTab$fdr<- p.adjust(outTab$pValue, method = "fdr")
outTabs<-outTab[outTab$pValue<0.05,]
TF_heatmap<-aucsmenas[outTabs$gene,] 
TF_heatmap<-`rownames<-`(apply(TF_heatmap,2,as.numeric), rownames(TF_heatmap))  

pdf(file.path(Figure3,"TF_scenic_fib.pdf"),height = 3.24,width=3.52)
pheatmap(TF_heatmap,
         scale = "row",
         show_colnames =T,
         show_rownames = T,
         cluster_cols = F,
         fontsize = 10,
         angle_col = c('45'),
         color=colorRampPalette(c('#1B9E77',"white",'#E41A1C'))(20)
)
dev.off()
#c('#984EA3','#377EB8','#1B9E77','#E41A1C')

