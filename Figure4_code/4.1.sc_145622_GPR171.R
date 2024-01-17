###scp for kegg of GPR171
library(dbplyr)
library(SCP)
library(Seurat)
library(BiocParallel)
library(data.table)

sce_LM_GSE178318<-NULL

RNA.res<-readRDS("~/Projects/ImmuneProGPR/Data/RNA.res_Zhangzeming_CancerCell_LM.RDS")
sce<-RNA.res_Zhangzeming_CancerCell_LM$RNA
colnames(sce@meta.data)
table(sce$assign.ident)
sce$assign.ident<-ifelse(sce$assign.ident=="CD8Tex","CD8T",sce$assign.ident)
sce$assign.ident<-ifelse(sce$assign.ident=="CD4Tconv","CD4T",sce$assign.ident)
sce$assign.ident<-ifelse(sce$assign.ident=="Mono/Macro","Macro",sce$assign.ident)
sce<-sce[,sce$assign.ident!="Others"]
table(sce$assign.ident)
###
mycol <- c("#223D6C","#D20A13","#FFD121","#088247",
           "#11AA4D","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C",
           "#E0367A","#D8D155","#64495D","#7CC767")
VlnPlot(object = sce,features="GPR171",group.by = "assign.ident", pt.size = 0.1,cols=mycol) & NoLegend()
FeatureStatPlot(sce, stat.by = "GPR171", group.by =  "assign.ident",pt.size = 0.1,plot_type = "bar",add_trend = TRUE)

FeatureStatPlot(sce, stat.by = "GPR171", group.by =  "assign.ident",pt.size = 0.1)

ggsave("GSE164522_GPR171bar.pdf")
####################################
cells_T<-c( "CD4T",        "CD8T",       "NK",     'Treg' )# "Mono/Macro", 
subsce<-sce[,sce$assign.ident %in% cells_T]
####

subsce$group<-ifelse(subsce@assays$RNA@counts['GPR171',]>0,1,0)
subsce$group<-paste0(subsce$assign.ident,"_",subsce$group)
table(subsce$group)
#subsce<-RunTSNE(subsce)
subsce<- RunUMAP(subsce, dims = 1:20)

####################
colnames(subsce@meta.data)
CellDimPlot(
  srt = subsce, group.by =  "assign.ident",
  reduction = "UMAP", theme_use = "theme_blank"
)
CellDimPlot(
  srt = subsce, group.by = "assign.ident", stat.by = "group",
  reduction = "TSNE", theme_use = "theme_blank"
)


sce<- RunUMAP(sce, dims = 1:20)
FeatureDimPlot(
  srt =sce, features = c("GPR171"),  
  reduction = "UMAP", theme_use = "theme_blank"
)
ggsave("GSE164522_GPR171bar_features.pdf")
CellDimPlot(
  srt = sce, group.by =  "assign.ident",
  reduction = "UMAP", theme_use = "theme_blank"
)
ggsave("GSE164522_GPR171bar_celltypes.pdf")

register(MulticoreParam(workers =20, progressbar = TRUE))


Degs<-function(x){
  print(x)
subsubsce<-subsce[,subsce$assign.ident==x]
Idents(subsubsce)<-subsubsce$group
sce.markers_GPR <- FindAllMarkers(subsubsce, 
                                  only.pos = TRUE, 
                                  min.pct = 0.15, logfc.threshold = 0.25)
subDEGs <- sce.markers_GPR #%>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
return(subDEGs)
}
subDEGslist<-lapply(cells_T,Degs)
subDEGslist<-rbindlist(subDEGslist)
DEGs<-subDEGslist
DEGs<-DEGs[-1,]
#subsce<-subsce[,]
#subsce<- RunDEtest(srt =subsce, group_by = "group", fc.threshold = 1, only.pos = FALSE)
#VolcanoPlot(srt = subsce, group_by = "group")
#ggsave(file = file.path(Figure_jbn,"VolcanoPlot_gene.pdf"))
#DEGs <- subsce@tools$DEtest_group$AllMarkers_wilcox
#DEGs <- DEGs[with(DEGs, avg_log2FC >0.1 & p_val_adj < 0.00001), ]
# Annotate features with transcription factors and surface proteins
subsce <- AnnotateFeatures(subsce, species = "Homo_sapiens", db = c("TF", "CSPA"))
factorgroup<-data.frame(table(DEGs$cluster))
DEGs$cluster<-factor(DEGs$cluster,levels= c( #"CD4T_0","CD4T_1" ,
                                             #"Mono/Macro_0", "Mono/Macro_1",
                                             "NK_0"  , "NK_1" , 
                                             "CD8T_0" , "CD8T_1", 
                                             "Treg_0","Treg_1"  ))
subsce$group<-factor(subsce$group,levels =c( "CD4T_0","CD4T_1" ,    
                                             #"Mono/Macro_0", "Mono/Macro_1",
                                             "NK_0"  , "NK_1" , 
                                             "CD8T_0" , "CD8T_1",
                                             "Treg_0","Treg_1"  ))
table(subsce$group)
##

mycol <-list(c("black","#D20A13","black","#D20A13","black","#D20A13","black","#D20A13","black","#D20A13","black","#D20A13"))

DefaultAssay(subsce)<-"RNA"
ht <- FeatureHeatmap(
  srt = subsce, 
  group.by = "group", 
  group_palette = "Paired",
  group_palcolor = mycol ,
  features = DEGs$gene, 
  feature_split = DEGs$cluster,
  feature_split_palcolor = mycol ,
  feature_annotation_palette = "Dark2",
  padjustCutoff = 0.05,
  #feature_annotation_palcolor =  mycol,
  species = "Homo_sapiens", 
  #db = c("GO_BP"，"KEGG"), #,
  db =  "GO_BP",
  anno_keys = TRUE, anno_features = TRUE,
  anno_terms = TRUE,
  use_raster=T,
  nlabel = 0,
  feature_annotation = c("TF", "CSPA"), 
  feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  #height = 9, 
  width = 3
)
print(ht$plot)


ggsave(ht$plot,file = "Figure3/function_KEGG_gene_heatmap.pdf")

table(subsce$assign.ident)

subsubsce<-subsce[,subsce$assign.ident=="CD8T"]
Idents(subsubsce)<-subsubsce$group
subsubsce$group<-factor(subsubsce$group)
table(subsubsce$group)
subsubsce<- RunDEtest(srt =subsubsce, group_by = "group", fc.threshold = 1, only.pos = FALSE)

subsubsce <- RunGSEA(
  srt = subsubsce, group_by = "group", db = c("GO_BP","KEGG", "WikiPathway", "Reactome"), species = "Homo_sapiens",
  DE_threshold = "p_val_adj < 0.05",
  simplify_similarityCutoff = 0.05
)
#c("KEGG", "WikiPathway", "Reactome", "PFAM", "MP"),
#results<-subsubsce@tools$GSEA_group_wilcox$results$`CD8T_1-GO_BP`@result

#results<-subsubsce@tools$GSEA_group_wilcox$results$`CD8T_1-WikiPathway`@result
#results<-results[results$Description %like% "immune",]
#results<-results[results$p.adjust<0.05,]
GSEAPlot(srt = subsubsce, group_by = "group", group_use = "CD8T_0")

GSEAPlot(subsubsce, db = "GO_BP", group_by = "group",  plot_type = "comparison")
GSEAPlot(subsubsce, db = "GO_BP", group_by = "group", group_use = "CD8T_1",  plot_type = "line")




#de_df<-DEGs[DEGs$cluster %in% c("CD8T_0","CD8T_1"),]
gsea_out <- RunGSEA(geneID = DEGs$gene, geneScore =  DEGs$avg_log2FC, 
                    geneID_groups = DEGs$cluster, db = "Reactome", species = "Homo_sapiens")
GSEAPlot(res = gsea_out, db = "Reactome", plot_type = "comparison")

gsea_out <- RunGSEA(geneID = DEGs$gene, geneScore =  DEGs$avg_log2FC, 
                    geneID_groups = DEGs$cluster, db = "GO_BP", species = "Homo_sapiens")
GSEAPlot(res = gsea_out, db = "GO_BP", plot_type = "comparison")

gsea_out <- RunGSEA(geneID = DEGs$gene, geneScore =  DEGs$avg_log2FC, 
                    geneID_groups = DEGs$cluster, db =c("GO_BP","KEGG", "WikiPathway", "Reactome"), species = "Homo_sapiens")
GSEAPlot(res = gsea_out, db = c("Reactome"),plot_type = "comparison")
?GSEAPlot


results<-subsubsce@tools$Enrichment_group_wilcox$results$`CD8T_1-GO_BP`@result
#results<-results[results$Description %like% "immune",]
results<-results[results$p.adjust<0.001,]
write.csv(results,file="CD8+GPR171.GO.csv")

results<-subsubsce@tools$Enrichment_group_wilcox$results$`CD8T_0-GO_BP`@result
#results<-results[results$Description %like% "immune",]
results<-results[results$p.adjust<0.000001,]
write.csv(results,file="CD8_lowGPR171.GO.csv")

EnrichmentPlot(
  srt =  subsubsce, group_by = "group", group_use = "CD8T_1",
  plot_type = "enrichmap"
)
EnrichmentPlot(
  srt =  subsubsce, group_by = "group", group_use = "CD8T_0",
  plot_type = "enrichmap"
)

subsubsce <- RunEnrichment(
  srt =  subsubsce, group_by = "group", db = "GO_BP",
  #db = c("KEGG", "WikiPathway", "Reactome", "PFAM", "MP"),
  species ="Homo_sapiens",# "Mus_musculus",
  DE_threshold = "avg_log2FC > log2(1) & p_val_adj < 0.001"
)


EnrichmentPlot(
  srt =  subsce, group_by = "newGroup", db = "MSigDB", group_use = c("Hepatocyte_0", "Hepatocyte_1"),
  plot_type = "bar"
)


EnrichmentPlot(
  srt =  subsce, group_by = "newGroup",  db = "MSigDB",group_use = c("Hepatocyte_0", "Hepatocyte_1"),
  plot_type = "enrichmap"
)
###
subsce<- RunUMAP(subsce, dims = 1:10)
subsce<- RunTSNE(object = subsce, dims = 1:10, do.fast = TRUE)#

subsce <- RunSlingshot(srt = subsce, group.by = "newGroup", reduction = "TSNE")
FeatureDimPlot(subsce, features = paste0("Lineage", 1:3), reduction = "TSNE", theme_use = "theme_blank")

##


immunegenes<-read.delim2("/home/data/gaoyuzhen/Projects/ImmuneProGPR/Results/ImmunecomparisonGene2.txt")

subsce$group
library(ggplot2)
ht <- GroupHeatmap(
  srt =subsce[,subsce$assign.ident=="CD8T"],
  features =immunegenes$Genes[immunegenes$types=="Co-inhibitors"],
  #palcolor =c("black","#E6AB02", "red"),
  #group_palette =NULL,
  #group_palcolor =c("#E6AB02", "red"),
  group.by =  "group",
  heatmap_palette = "YlOrRd",
  #cell_annotation = c( "cohort"),
  #cell_annotation_palette = c( "Paired"),
  show_row_names = TRUE, row_names_side = "left",
  #add_violin = TRUE,
  add_bg = TRUE,
  nlabel = 0,
  add_dot = TRUE, add_reticle = TRUE,
  border=TRUE
)
print(ht$plot)
ggsave("Figure3/immunecheckponts_T_GPR171.pdf",width = 4.9,height = 9.93)


ht2 <- GroupHeatmap(
  srt =subsce[,subsce$assign.ident=="CD8T"],
  features =immunegenes$Genes[immunegenes$types!="Co-inhibitors"] ,
  #palcolor =c("black","#E6AB02", "red"),
  #group_palette =NULL,
  #group_palcolor =c("#E6AB02", "red"),
  group.by =  "group",
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
print(ht2$plot)

ht$plot+ht2$plot
ggsave("Figure2.1/immunecheckponts_T.pdf",width = 4.9,height = 9.93)