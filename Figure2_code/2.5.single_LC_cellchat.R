
library(data.table)
library(ggsci)
library(pheatmap)
library(AUCell)
library(GEOquery)
library(tidyverse)
library(Seurat)
library(RColorBrewer)
library(dplyr)
library(Seurat)
library(ggplot2)
library(umap)
library(hypeR)
library(msigdbr)
library(pheatmap)
library(corrplot)#


#load("LM_meta_sce.Rdata")
#subsce<-readRDS("LM_subsce.RDS")
#/home/data/gaoyuzhen/Projects/LiverCancer/singlecellData/LIHC_GSE140228_10X_expression.h5
#


data.input<-Read10X_h5("/home/data/gaoyuzhen/Projects/LiverCancer/singlecellData/LIHC_GSE166635_expression.h5", use.names = TRUE, unique.features = TRUE)
phe<- read.delim2("/home/data/gaoyuzhen/Projects/LiverCancer/singlecellData/LIHC_GSE166635_CellMetainfo_table.tsv")

table(phe$Celltype..major.lineage.)

meta<-phe

colnames(meta)
meta$id<-meta$Cell
meta<-meta[,c("id","Celltype..major.lineage.")]
colnames(meta)<-c("Index","Cell_type")
rownames(meta)<-meta$Index


table(meta$Cell_type)
meta<-meta[!meta$Cell_type %in% c("Endothelial","Epithelial", "Fibroblasts"),]
#################################
#save(meta,file="LM_meta_sce.Rdata")
#saveRDS(subsce,file="LM_subsce.RDS")
###
save(meta,file="LC_meta_sce.Rdata")###


library(data.table)
library(ggsci)
library(pheatmap)
library(AUCell)
library(GEOquery)
library(tidyverse)
library(Seurat)
library(RColorBrewer)
library(dplyr)
library(Seurat)
library(ggplot2)
library(umap)
library(hypeR)
library(msigdbr)
library(pheatmap)
library(corrplot)#
load("LC_meta_sce.Rdata")
data.input<-Read10X_h5("/home/data/gaoyuzhen/Projects/LiverCancer/singlecellData/LIHC_GSE166635_expression.h5", use.names = TRUE, unique.features = TRUE)
data.input<-data.input[,colnames(data.input) %in% meta$Index]
data.input<-as.matrix(data.input)

###
library(CellChat)
library(patchwork)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(SeuratData)
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "Cell_type")
###
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "Cell_type") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
###我们的数据库 CellChatDB 是一个手动整理的文献支持的配体受体在人和小鼠中的交互数据库。
##小鼠中的CellChatDB包含2，021个经验证的分子相互作用，包括60%的自分泌/旁分泌信号相互作用
##、21%的细胞外基质（ECM）受体相互作用和19%的细胞-细胞接触相互作用。人的CellChatDB包含1，939个经验证的分子相互作用，
##包括61.8%的自分泌/旁分泌信号相互作用、21.7%的细胞外基质（ECM）受体相互作用和16.5%的细胞-细胞接触相互作用。
##
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
#showDatabaseCategory(CellChatDB)
#### Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use
##预处理
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 70) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
future::plan("multiprocess", workers = 70) # do parallel
cellchat <- identifyOverExpressedInteractions(cellchat)
future::plan("multiprocess", workers = 70) # do parallel
cellchat <- projectData(cellchat, PPI.human)
##在分析未分类的单细胞转录组时，假设丰富的细胞群倾向于发送比稀有细胞群更强的信号，
#CellChat 还可以在概率计算中考虑每个细胞组中细胞比例的影响。用户可以设置population.size = TRUE
#计算通信概率并推断cellchat网络##type = "truncatedMean"和对trim = 0.1。
future::plan("multiprocess", workers = 70) # do parallel
cellchat <- computeCommunProb(cellchat, raw.use = TRUE,population.size = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
future::plan("multiprocess", workers = 70) # do parallel
cellchat <- filterCommunication(cellchat, min.cells = 2)
saveRDS(cellchat,file="cellchat.rds")##first 

#####################

cellchat<-readRDS("cellchat.rds")
##
##
df.net <- subsetCommunication(cellchat)
##
write.csv(df.net,file="df.net_all_LC.csv")
#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))#将推断的细胞-细胞通信从细胞组1和2发送到细胞组4和5。
#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))#通过向WNT和TGFb发出信号来调节推断的细胞通信。
##在信号通路级别推断细胞-细胞通信
cellchat <- computeCommunProbPathway(cellchat)
##计算整合的细胞通信网络 我们可以通过计算链接数或汇总通信概率来计算整合的细胞通信网络。用户还可以通过设置sources.use和targets.use`
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count,
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= T,
                 title.name = "Number of interactions")

ggsave("full_manuscript.pdf")
netVisual_circle(cellchat@net$weight,
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength")
##
colnames(mat)
mat <- cellchat@net$count
mat<-mat[!colnames(mat) %in% c("Endothelial","Epithelial", "Fibroblasts"),!rownames(mat) %in% c("Endothelial","Epithelial", "Fibroblasts")]
####

groupSize <- as.numeric(table(cellchat@idents[!cellchat@idents %in% c("Endothelial","Epithelial", "Fibroblasts")]))

groupSize<-groupSize[groupSize!=0]
###
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:nrow(mat)) { 
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat)) 
  mat2[, i] <- mat[,i]  
  pdf(paste0(gsub("/",".",rownames(mat))[i],"LiverCancer_network.pdf"))
  
  netVisual_circle(mat2, vertex.weight = groupSize, label.edge= T,  weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  dev.off()
}

###########ggplot
Cellchat_results_comparison <- read_excel("Cellchat_results_comparison.xlsx",1)
Cellchat_results_comparison <- read_excel("Cellchat_results_comparison.xlsx",2)
data<-Cellchat_results_comparison
colnames(data)
data$x<-data$...1

library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(ggpubr)
#install.packages("hrbrthemes", repos = "https://cinc.rud.is") 


ggplot(data) +
  geom_segment( aes(x=x, xend=x, y=PLC, yend=LM), color="grey40") +
  geom_point( aes(x=x, y=PLC), color="#008280CC", size=3 ) +
  geom_point( aes(x=x, y=LM), color="#ED0000CC", size=3 ) +
  coord_flip()+
  #theme_ipsum() +
  geom_hline(aes(yintercept=8),color="black")+
  theme(
    legend.position = "none",
  ) +
  xlab("") +
  ylab("The CellPhone in immunecells with Malignant")
ggsave("figure2_cellphone1.pdf")

ggsave("figure2_cellphone2.pdf")

