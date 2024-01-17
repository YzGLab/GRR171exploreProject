##泛癌相关性 分析！！！ 类似于 GPR171 基因
###相关性
load("~/Projects/MET500Projects/immunefeaturs.Rdata")
load("~/Projects/MET500Projects/phe_met.Rdata")
load("~/Projects/MET500Projects/metExp500.Rdata")

#########################

#devtools::install_local("/home/data/gaoyuzhen/Projects/ImmunePro/ImmuCellAI-main.zip")
library("ImmuCellAI")
phe_met$group<-factor(phe_met$group)
rownames(phe_met)<-phe_met$run.id
ai<-ImmuCellAI_new(metExp,data_type="rnaseq",group_tag = phe_met$group,response_tag=phe_met$group)
#GPR<-metExp["GPR171",] %>% t() %>% cbind(Sample_id=rownames(.),.) %>% data.frame()
GPR<-metExp[c("GPR171","CD8A"),] %>% t() %>% cbind(Sample_id=rownames(.),.) %>% data.frame()
GPR$GPR171<-as.numeric(GPR$GPR171)/as.numeric(GPR$CD8A)
#rt<-merge(GPR,rt,by="Sample_id")
#GPR$Sample_id<-gsub("-",".",GPR$Sample_id)
#计算不同人群的相关性
library(data.table)
library(ggplot2)
library(limma)
library(scales)
library(ggtext)
library(reshape2)
library(tidyverse)
library(ggpubr)
#devtools::install_github("Hy4m/linkET", force = TRUE)
library(ggprism)
library(forcats)
library(ggh4x)
library(ggplot2)
library(ggplot2)
library(ggfun)
library(splines)
library(linkET)

KR <- as.data.frame(immPath.score[,1])
colnames(KR) <- "GPR171"
immPath.score2 <- immPath.score[-1,]
immCorsigene <- NULL
for (i in colnames(immPath.score)) {
  cr <- cor.test(as.numeric(immPath.score[,i]),
                 as.numeric(KR[,1]),      #选择某一行作为link的点，按需更改 
                 method = "pearson")
  immCorsigene <- rbind.data.frame(immCorsigene,
                                   data.frame(gene = "GPR171",     #link点的行名，按需更改
                                              path = i,
                                              r = cr$estimate,
                                              p = cr$p.value,
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
}
immCorsigene$sign <- ifelse(immCorsigene$r > 0,"pos","neg")
immCorsigene$absR <- abs(immCorsigene$r)
immCorsigene$rSeg <- as.character(cut(immCorsigene$absR,c(0,0.25,0.5,0.75,1),labels = c("0.25","0.50","0.75","1.00"),include.lowest = T))
immCorsigene$pSeg <- as.character(cut(immCorsigene$p,c(0,0.001,0.01,0.05,1),labels = c("<0.001","<0.01","<0.05","ns"),include.lowest = T))
immCorsigene[nrow(immCorsigene),"pSeg"] <- "Not Applicable"

immCorsigene$rSeg <- factor(immCorsigene$rSeg, levels = c("0.25","0.50","0.75","1.00"))
immCorsigene$pSeg <- factor(immCorsigene$pSeg, levels = c("<0.001","<0.01","<0.05","Not Applicable","ns"))
immCorsigene$sign <- factor(immCorsigene$sign, levels = c("pos","neg"))
colnames(immCorsigene)

#immCorsigene<-immCorsigene[,colnames(immCorsigene) %in% immCorsigene$path]

W<-apply(immPath.score,2,as.numeric)
W<-as.matrix(W)
library(corrplot)
p4 <- qcorrplot(correlate(W), 
                type = "lower",
                diag = F) + 
  geom_square() +
  geom_couple(data = immCorsigene, 
              mapping = aes(colour = pSeg, size = rSeg, linetype = sign,x = .x - 4, ),
              curvature = nice_curvature(0.15),
              node.colour = c("blue", "#F3C3B6"),
              node.fill = c("#78c8e2", "#35825a"),
              node.size = c(5, 3)) +
  scale_size_manual(values = c(0.5, 1, 1.5, 2)) +
  geom_mark(# 添加r值与显著性标�?
    sep = '\n', 
    size = 3, 
    sig_level = c(0.05, 0.01, 0.001), # 显著性水平设�?
    sig_thres = 0.001 # 显著性阈值，p值大于阈值的相关性系数不会被绘制�?
  )+
  theme(legend.background=element_roundrect(color="#808080", linetype=2))+
  # geom_number(data = get_data(p.value < 0.05, type = "lower"),aes(num = r),size = 2,)+   #显示相关系数R=�?
  scale_color_manual(values = c("#006D77","#83C5BE","#669BBC","#FFDDD2","#E29578")) +
  scale_fill_gradient2(low = "black",mid = "white",high = "#E31A1C",midpoint=0)+
  guides(colour = guide_legend(title = "pSeg", 
                               override.aes = list(size = 3), 
                               order = 1))


p4


ggsave(file="Figure3/Tcell_cor_GPR_LM.pdf",width = 12.8,height = 10)
##################################