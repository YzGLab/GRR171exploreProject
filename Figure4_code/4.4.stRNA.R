############
##########
library(SPATA2)
library(SPATAData)
library(tidyverse)

###
LM1_spatial <-
  initiateSpataObject_10X(
    directory_10X = '/home/data/gaoyuzhen/Projects/ImmuneProGPR/Data/CancerDiscovery_livermeta/ST/ST_liver1', # the directory from which to load the data
    sample_name = "ST_liver1"
  )
getDirectoryInstructions(object = LM1_spatial)
# set/change the current default directory
LM1_spatial <- setSpataDir(LM1_spatial, dir = "Data/spata_objects/LM1_spatial.RDS")
# quickly save the `spata2` object under the default directory 
saveSpataObject(object = LM1_spatial)
##############################
LM1_spatial<-readRDS("Data/spata_objects/LM1_spatial.RDS")
LM1_spatial <- setActiveMatrix(object = LM1_spatial, mtr_name = "scaled")
##机器学习去除噪声 
LM1_spatial <-
  runAutoencoderDenoising(
    object = LM1_spatial, 
    activation = "selu", 
    bottleneck = 56, 
    epochs = 20, 
    #layers = c(128, 64, 32), 
    dropout = 0.1
  )

getExpressionMatrixNames(object =LM1_spatial)
# active expression matrix after denoising
getActiveMatrixName(object = LM1_spatial)
###############################################################
# discard the one added above, to create new one
#LM1_spatial <- discardFeatures(object = LM1_spatial, feature_names = "histology")
##################################
LM1_spatial <- createSpatialSegmentation(LM1_spatial)
#################################################
LM1_spatial@fdata$ST_liver1$histology<-ifelse(LM1_spatial@fdata$ST_liver1$histology=="tumor","tumor","non-tumor")
#install.packages("Cairo",force=T)
plotSurface(object =LM1_spatial, color_by = "histology", pt_clrp= "turbo",clrp_adjust = c("non-tumor" ="#00a2b3", "tumor" =  "#cf3e53"))

ggsave(file = file.path(Figure3,"C1_plotCnvHeatmap_histology.pdf"),width =11.97,height = 6.5)####
# names of grouping variables
getGroupingOptions(object = LM1_spatial)
# only histology
plotSurface(object =LM1_spatial, pt_alpha = 0,bcsp_rm=T)
###############################展示单基因
LM1_spatial <- setActiveMatrix(object = LM1_spatial, mtr_name = "denoised")
###############################展示多基因
plotSurfaceComparison(
  object = LM1_spatial, 
  #color_by = c("CD274", "CD8A","GPR171","CD68","PDCD1","FOXP3","SPP1","CD4"), 
  color_by = c("CD274", "CD8A","GPR171","PDCD1","FOXP3","CD4"), 
  pt_clrsp = "Greens 3", 
  display_image = TRUE, 
  smooth = TRUE, 
  alpha_by = TRUE
) 

ggsave(file = "/home/data/gaoyuzhen/Projects/ImmuneProGPR/Figure3/LM1_SurfaceComparison_genes.pdf",width =15.97,height =11.5)####
##############################
getGroupingOptions(LM1_spatial)
plotSurface(
  object = LM1_spatial, 
  color_by = "seurat_clusters"
  #pt_clrp = "uc"
)
###
LM4_spatial <-
  initiateSpataObject_10X(
    directory_10X = '/home/data/gaoyuzhen/Projects/ImmuneProGPR/Data/CancerDiscovery_livermeta/ST/ST_liver4', # the directory from which to load the data
    sample_name = "ST_liver4"
  )
getDirectoryInstructions(object = LM4_spatial)
# set/change the current default directory
LM4_spatial <- setSpataDir(LM4_spatial, dir = "Data/spata_objects/LM4_spatial.RDS")
# quickly save the `spata2` object under the default directory 
saveSpataObject(object = LM4_spatial)
##############################
LM4_spatial <- setActiveMatrix(object = LM4_spatial, mtr_name = "scaled")
##机器学习去除噪声 
LM4_spatial <-
  runAutoencoderDenoising(
    object = LM4_spatial, 
    activation = "selu", 
    bottleneck = 56, 
    epochs = 20, 
    #layers = c(128, 64, 32), 
    dropout = 0.1
  )

getExpressionMatrixNames(object =LM4_spatial)
# active expression matrix after denoising
getActiveMatrixName(object = LM4_spatial)
###############################################################
# discard the one added above, to create new one
#LM4_spatial <- discardFeatures(object = LM4_spatial, feature_names = "histology")
##################################
LM4_spatial <- createSpatialSegmentation(LM4_spatial)
#################################################
LM4_spatial@fdata$ST_liver4$histology<-ifelse(LM4_spatial@fdata$ST_liver4$histology=="tumor","tumor","non-tumor")
#install.packages("Cairo",force=T)
plotSurface(object =LM4_spatial, color_by = "histology", pt_clrp= "turbo",clrp_adjust = c("non-tumor" ="#00a2b3", "tumor" =  "#cf3e53"))

ggsave(file = file.path(Figure3,"C1_plotCnvHeatmap_histology.pdf"),width =11.97,height = 6.5)####
# names of grouping variables
getGroupingOptions(object = LM4_spatial)
# only histology
plotSurface(object =LM4_spatial, pt_alpha = 0,bcsp_rm=T)
###############################展示单基因
LM4_spatial <- setActiveMatrix(object = LM4_spatial, mtr_name = "denoised")
###############################展示多基因
plotSurfaceComparison(
  object = LM4_spatial, 
  #color_by = c("CD274", "CD8A","GPR171","PDCD1","CD4"), 
  color_by = c("CD274", "CD8A","GPR171","PDCD1","FOXP3","CD4"), 
  pt_clrsp = "Greens 3", 
  display_image = TRUE, 
  smooth = TRUE, 
  alpha_by = TRUE
) 

ggsave(file = "/home/data/gaoyuzhen/Projects/ImmuneProGPR/Figure3/LM4_SurfaceComparison_genes.pdf",width =15.97,height =11.5)####


###################

###
C1_spatial <-
  initiateSpataObject_10X(
    directory_10X = '/home/data/gaoyuzhen/Projects/ImmuneProGPR/Data/CancerDiscovery_livermeta/ST/ST-colon1', # the directory from which to load the data
    sample_name = "ST_colon1"
  )
getDirectoryInstructions(object = C1_spatial)
# set/change the current default directory
C1_spatial <- setSpataDir(C1_spatial, dir = "Data/spata_objects/C1_spatial.RDS")
# quickly save the `spata2` object under the default directory 
saveSpataObject(object = C1_spatial)
##############################
C1_spatial <- setActiveMatrix(object = C1_spatial, mtr_name = "scaled")
##机器学习去除噪声 
C1_spatial <-
  runAutoencoderDenoising(
    object = C1_spatial, 
    activation = "selu", 
    bottleneck = 56, 
    epochs = 20, 
    #layers = c(128, 64, 32), 
    dropout = 0.1
  )

getExpressionMatrixNames(object =C1_spatial)
# active expression matrix after denoising
getActiveMatrixName(object = C1_spatial)
###############################################################
# discard the one added above, to create new one
#C1_spatial <- discardFeatures(object = C1_spatial, feature_names = "histology")
##################################
C1_spatial <- createSpatialSegmentation(C1_spatial)
#################################################
C1_spatial@fdata$ST_colon1$histology<-ifelse(C1_spatial@fdata$ST_colon1$histology=="tumor","tumor","non-tumor")
#install.packages("Cairo",force=T)
plotSurface(object =C1_spatial, color_by = "histology", pt_clrp= "turbo",clrp_adjust = c("non-tumor" ="#00a2b3", "tumor" =  "#cf3e53"))

ggsave(file = file.path(Figure3,"C1_plotCnvHeatmap_histology.pdf"),width =11.97,height = 6.5)####
# names of grouping variables
getGroupingOptions(object = C1_spatial)
# only histology
plotSurface(object =C1_spatial, pt_alpha = 0,bcsp_rm=T)
###############################展示单基因
C1_spatial <- setActiveMatrix(object = C1_spatial, mtr_name = "denoised")
###############################展示多基因
plotSurfaceComparison(
  object = C1_spatial, 
  color_by = c("CD274", "CD8A","GPR171"), 
  pt_clrsp = "Greens 3", 
  display_image = TRUE, 
  smooth = TRUE, 
  alpha_by = TRUE
) 

ggsave(file = "LM1_SurfaceComparison_genes.pdf",width =8.97,height = 6.5)####