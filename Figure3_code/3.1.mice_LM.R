#####mouse marker 
library(readxl)
library(ggsci);library(pheatmap)
library(Seurat)
library(RColorBrewer);library(corrplot)#
library(reshape2)
library(ggplot2)
library(data.table);library(tidyr)
library(dplyr)
library(clustree);
library(gplots);library(ComplexHeatmap)
library(ggpubr);library(scales)
library(RColorBrewer);library(MAESTRO)
library(tidyverse)
library(future)

mouseSce<-readRDS("/home/data/gaoyuzhen/Projects/Forxia_macDICER/sc_liver.RDS")
library(scRNAtoolVis)
library(jjPlot)
mouseSce<-mouseSce[,mouseSce$celltype!=" "]
#pdf(file.path(Figure4,"metabolism_DimPlot_all_site_mac_total.pdf"),width = 5,height =4)
scRNAtoolVis::scatterCellPlot(object =mouseSce,rm.axis = F,cell.id = "celltype", dim  = "umap",color = mycol[c(1:8)] )
scRNAtoolVis::scatterCellPlot(object =mouseSce)

library(grid)
library(ggplot2)
library(SCP)       
library(Seurat)
table(mouseSce$orig.ident)
mouseSce$celltype[mouseSce$celltype=="Natural killer cell"]="T cell"
mouseSce$source<-ifelse(mouseSce$orig.ident=="GSM4771486","sc","sc_liver")
#mouseSce<-Standard_SCP(mouseSce)


CellDimPlot(
  srt = mouseSce, group.by =  "celltype",
  palcolor =c('#1B9E77','#377EB8','#984EA3','#E41A1C','#FFFFB3','#FDBF6F',"#7A142C"),
  split.by = "source",
  cells.highlight = TRUE,
  reduction = "UMAP", theme_use = "theme_blank"
)

ggsave("Figure2.1/CellDimPlot_mouse.pdf")



CellStatPlot( mouseSce,stat.by ="celltype", group.by ="source",stat_type = "percent",position = "dodge",
              palcolor =c('#1B9E77','#377EB8','#984EA3','#E41A1C','#FFFFB3','#FDBF6F',"#7A142C"), label = TRUE)
ggsave("Figure2.1/CellDimPlot_mouse_ratio_celltype.pdf")
#CellStatPlot((mouseSce,stat.by ="source", group.by ="celltype",label = TRUE) %>% panel_fix(height = 2, width = 3)
CellStatPlot(mouseSce,stat.by ="source", group.by ="celltype", stat_type = "percent", 
             palcolor = c("black",'#E41A1C'), 
             position = "dodge", label = TRUE)
ggsave("Figure2.1/CellDimPlot_mouse_ration.pdf")
CellStatPlot(srt = mouseSce, stat.by = "source",  group.by ="celltype",palcolor = c("black",'#E41A1C'),  label = TRUE)
ggsave("Figure2.1/CellDimPlot_mouse_ratio2.pdf")


CellStatPlot( mouseSce,stat.by ="celltype", group.by ="source",plot_type = "area")
CellStatPlot(mouseSce,stat.by ="celltype",palcolor =c('#1B9E77','#377EB8','#984EA3','#E41A1C','#FFFFB3','#FDBF6F',"#7A142C"),  group.by ="source", plot_type = "dot")

FeatureDimPlot(mouseSce,features = "Gpr171",split.by ="source")

Gpr171<-mouseSce@assays$RNA@counts["Gpr171",]
Gpr171<-ifelse(Gpr171>0,1,0)
mouseSce<-AddMetaData(mouseSce,Gpr171,col.name = "Gpr171")
FeatureStatPlot(mouseSce, stat.by = "Gpr171", group.by = "celltype",split.by = "source")
VlnPlot(mouseSce,features = "Gpr171",group.by = "celltype",split.by = "source" )
###############
table(mouseSce$celltype)
DefaultAssay(mouseSce)<-"RNA"
Idents(mouseSce)<-mouseSce$source
sce.markers2<-FindAllMarkers(mouseSce[,mouseSce$celltype=="T cell"],
                             only.pos = TRUE, 
                             min.pct = 0.05,
                             min.diff.pct = 0.05)
sce.markers2["Gpr171",]
Tsce<-mouseSce[,mouseSce$celltype=="T cell"]
Tsce <- RunDEtest(srt =Tsce, group_by = "source", fc.threshold = 1, only.pos = FALSE)
table(Tsce$celltype,Tsce$source)
VolcanoPlot(srt =Tsce, group_by ="source")
ggsave("Figure2.1/CellDimPlot_mouse_VolcanoPlot_T.pdf")


VlnPlot(Tsce,features = "Gpr171",group.by="source" )
ggsave("Figure2.1/VlnPlot_Gpr171_T.pdf")


DefaultAssay(Tsce)<-"integrated"
Tsce <- FindNeighbors(Tsce, dims = 1:15)
Tsce <- FindClusters(Tsce, resolution = 0.1)
set.seed(123)
Tsce<-RunUMAP(Tsce,dims = 1:15, do.fast = TRUE)
Tsce <- RunTSNE(object = Tsce, dims = 1:15, do.fast = TRUE)
colnames(Tsce@meta.data)

CellDimPlot(
  srt = Tsce, 
  group.by =  "seurat_clusters",
  palcolor =c('#1B9E77','#377EB8','#984EA3','#E41A1C','#FFFFB3','#FDBF6F',"#7A142C"),
  #split.by = "source",
  reduction = "UMAP", theme_use = "theme_blank"
)

Tsce$group<-ifelse(Tsce$seurat_clusters==0,1,0)

CellStatPlot(Tsce,stat.by ="group",palcolor =c('#1B9E77','#377EB8','#984EA3','#E41A1C','#FFFFB3','#FDBF6F',"#7A142C"),  group.by ="source", plot_type = "dot")

CellStatPlot( Tsce,stat.by ="group", group.by ="source",stat_type = "percent",position = "dodge",
              palcolor =c('#1B9E77','#377EB8','#984EA3','#E41A1C','#FFFFB3','#FDBF6F',"#7A142C"), label = TRUE)

CellStatPlot(Tsce,stat.by ="source", group.by ="seurat_clusters", stat_type = "percent", 
             palcolor = c("black",'#E41A1C'), 
             position = "dodge", label = TRUE)

DefaultAssay(Tsce)<-"RNA"
subTsce <- RunDEtest(srt =Tsce,group_by = "source", fc.threshold = 1, only.pos = FALSE)
Idents(subTsce)<-subTsce$source
subTsce <- RunEnrichment(
  srt = subTsce, group_by =  "source", db = "KEGG", species = "Mus_musculus",
  DE_threshold = "p_val_adj < 0.05"
)

EnrichmentPlot(
  srt =subTsce, group_by = "source", db ="KEGG", group_use = c("sc", "sc_liver"),
  plot_type = "bar", topTerm = 10)

EnrichmentPlot(
  srt = subTsce, group_by = "source", db = "KEGG",  group_use = c("sc", "sc_liver"),
  plot_type = "network", topTerm = 10)
###


Tsce <- RunEnrichment(
  srt = Tsce, group_by ="source",
  db = c("KEGG", "WikiPathway", "Reactome"),
  db_combine = TRUE,
  species = "Mus_musculus"
)
EnrichmentPlot(Tsce, db = "Combined", group_by = "source", plot_type = "comparison")
################


DefaultAssay(Tsce)<-"RNA"
Idents(Tsce)<-Tsce$source
sce.markers3<-FindAllMarkers(Tsce,
                             only.pos = TRUE, 
                             min.pct = 0.25,
                             min.diff.pct = 0.25)
sce.markers3<-sce.markers3[sce.markers3$cluster!="sc",]
write.csv(sce.markers3,file='/home/data/gaoyuzhen/Projects/ImmuneProGPR/Figure2.1/sce.markers3_sc_liver_T.csv')

##pathway analysis
library(KEGGREST)  
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)


Myenrich_mus <- function(genes, category = c("kegg", "gobp"), 
                     geneid = c("SYMBOL", "ENTREZID", "ENSEMBL", "UNIPROT")){
  library(clusterProfiler)
  category <- match.arg(category)
  geneid <- match.arg(geneid)
  if (geneid != "ENTREZID"){
    genes <- bitr(genes, fromType = geneid,toType ="ENTREZID",OrgDb="org.Mm.eg.db") %>% .[, 2] %>% as.character()
  }else{genes <- genes}
  if (category == "kegg"){
    enrich <- enrichKEGG(genes,keyType = "kegg", organism = 'mmu',
                         pAdjustMethod = "BH",pvalueCutoff  = 0.001,
                         qvalueCutoff  = 0.001, maxGSSize = 5000)
  }
  if (category == "gobp"){
    enrich <- enrichGO(genes, OrgDb="org.Mm.eg.db",ont= "BP",pAdjustMethod = "BH",
                       pvalueCutoff  = 0.0001,qvalueCutoff  = 0.001, readable = TRUE)
  }
  return(enrich)
}

###########################
ProSig <- list(sce.markers3$gene)####设置基因集
Sigbp <- lapply(ProSig, Myenrich_mus , category = "gobp", geneid = "SYMBOL")
Sigkegg <- lapply(ProSig, Myenrich_mus , category = "kegg", geneid = "SYMBOL")


edox1 <- pairwise_termsim(Sigbp[[1]])
edox2 <- pairwise_termsim(Sigkegg[[1]])

p1 <- treeplot(edox1, hclust_method = "average")
p2 <- treeplot(edox2, hclust_method = "average")
p1+p2
ggsave(file.path(Figure2,"CRC_Tumor_vs_LM_Tumor_pathway_GO.pdf"),
       width = 11.5, height =4.73,
       bg = "transparent", dpi = 300, useDingbats=FALSE)

p1<-emapplot(edox1, hclust_method = "average")
p2<-emapplot(edox2, hclust_method = "average")
p1+p2
ggsave(file.path(Figure2,"CRC_Tumor_vs_LM_Tumor_pathway_GO_network.pdf"),
       width = 12.5, height =8.73,
       bg = "transparent", dpi = 300, useDingbats=FALSE)
emapplot(edox1, layout="kk")
emapplot(edox2, cex_category=1)
####################
library(biomaRt)
# 创建一个人类和mouse的 mart 对象
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
# 人类基因列表
human_genes <- immunegenes$Genes
# 使用 biomaRt 进行基因名称的转换
genes_converted <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol",
                          values = human_genes, mart = human,
                          attributesL = c("mgi_symbol"), martL = mouse)
# 输出转换后的小鼠基因名称
mouse_genes <- genes_converted[, 2]
print(mouse_genes)
genes_converted$Genes<-genes_converted$HGNC.symbol

genes_converted<-merge(genes_converted,immunegenes,by="Genes")
immunegenes<-genes_converted

mouse_genes_co<-mouse_genes[mouse_genes %in% immunegenes$Genes[immunegenes$types!="Co-inhibitors"] ]
immunegenes<-read.delim2("/home/data/gaoyuzhen/Projects/ImmuneProGPR/Results/ImmunecomparisonGene2.txt")

ht <- GroupHeatmap(
  srt = Tsce,
  features =immunegenes$MGI.symbol[immunegenes$types=="Co-inhibitors"],
  #palcolor =c("black","#E6AB02", "red"),
  #group_palette =NULL,
  #group_palcolor =c("#E6AB02", "red"),
  group.by =  "source",
  heatmap_palette = "YlOrRd",
  #cell_annotation = c( "cohort"),
  cell_annotation_palette = c( "Paired"),
  show_row_names = TRUE, row_names_side = "left",
  #add_violin = TRUE,
  add_bg = TRUE,
  nlabel = 0,
  add_dot = TRUE, add_reticle = TRUE,
  border=TRUE
)
print(ht$plot)



ggsave("Figure2.1/immunecheckponts_T.pdf")
