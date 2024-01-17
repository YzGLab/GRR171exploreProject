####1
library(tableone)
library(data.table)
library(tidyverse)
library(readxl)
library(readr)
library(survminer)
library(survival)
library(ggplot2)
library(tibble)
library(scales)
library(ggrepel)
library(forcats)
library(reshape2)
library(cowplot)
library(data.table)
library(tidyr)
####

load("livemeta_cohort.Rdata")
##data input
data_immunotherapy<- read_csv("/home/data/gaoyuzhen/Projects/ImmunePro/Figure1/mskcc_combined_updated.csv")
####################

TCGA<-read_excel("/home/data/gaoyuzhen/Projects/ImmunePro/Figure1/PANCANCER_information.xlsx")
total <- read_excel("/home/data/gaoyuzhen/Projects/ImmunePro/Figure1/CellData_1.xlsx")
colnames(total )
##
total<-total[total$sample_type=="Metastasis",]
data_immunotherapy$sample_id<-data_immunotherapy$ID
data_immunotherapy<-data_immunotherapy[data_immunotherapy$ID %in% total$sample_id,]
data_immunotherapy<-merge(data_immunotherapy,total,by="sample_id")
table(data_immunotherapy$MetaLiver)
################

tmp<-data_immunotherapy
tmp$MetaLiver<-ifelse(tmp$MetaLiver==3,"Other metastasis","liver metastasis")
shape_level <- nlevels(factor(tmp$Cancer_Type))
if (shape_level < 15){
  shapes = (0:shape_level) %% 15
} else{
  shapes = c(0:14,c((15:shape_level) %% 110 + 18))
}

colnames(tmp)
mycol<-c("#E31A1C","black","#33A02C","#BC80BD","#377EB8", "#6A3D9A" )

  tmp[,"tmb"]<-as.numeric(tmp[,"tmb"])
  topvealues<-median(tmp[,"tmb"])+3*sd(tmp[,"tmb"])
  tmp[,"tmb"]<-ifelse(tmp[,"tmb"]>  topvealues,  topvealues,tmp[,"tmb"])
  tmp$MetaLiver<-ifelse(tmp$MetaLiver=="liver metastasis",1,0)
  tmp$MetaLiver<-as.factor(tmp$MetaLiver)


  
  ################
  for(type in unique(tmp$Cancer_Type)){
    tmp1<-tmp[tmp$Cancer_Type==type,]
  c<- ggplot(tmp1,
             aes_string(x="MetaLiver", y="tmb", 
                        fill ="MetaLiver", 
                        color ="MetaLiver")) +
    geom_point(aes(shape=Cancer_Type), 
               position = position_jitter(width = .15,height=-0.7),
               size=1)+
    xlab("")+
    scale_shape_manual(values=shapes) +
    scale_color_manual(values=c(mycol)) +
    scale_fill_manual(values=c(mycol)) +
    geom_boxplot(notch = F, alpha = 0.65, width=0.2, outlier.size = 0.65)+
    geom_violin(position = position_nudge(x = 0),width=0.9,alpha = 0.4)+
    #ggtitle("CPSS score", "stratified by InternStudent") + 
    theme_light() +
    theme(axis.text.x = element_text(angle =30, hjust = 1,size = 10), 
          axis.text.y = element_text(colour="grey",angle = 90, size = 12),
          axis.title.y = element_text(angle = 90, size = 15)) +
    theme(legend.position = "right")+
    theme(axis.title=element_text(size=20),
          axis.text = element_text(face = "bold",size = 16),
          axis.ticks.length=unit(.4, "cm"),
          axis.ticks = element_line(colour = "black", size = 1),
          #panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1.5),
          plot.margin = margin(1, 1, 1, 1, "cm"))
  
  c+stat_compare_means(method = 'wilcox.test')
  ggsave(paste0("Figure1/TMB/",type,"TMB_LM.pdf"),height=6.5,width=6.5)
  }###
  
  
  
  c<- ggplot(tmp,
             aes_string(x="MetaLiver", y="tmb", 
                        fill ="MetaLiver", 
                        color ="MetaLiver")) +
    geom_point(aes(shape=Cancer_Type), 
               position = position_jitter(width = .15,height=-0.7),
               size=1)+
    xlab("")+
    scale_shape_manual(values=shapes) +
    scale_color_manual(values=c("black","#E31A1C")) +
    scale_fill_manual(values=c("black","#E31A1C")) +
    geom_boxplot(notch = F, alpha = 0.65, width=0.2, outlier.size = 0.65)+
    geom_violin(position = position_nudge(x = 0),width=0.9,alpha = 0.4)+
    #ggtitle("CPSS score", "stratified by InternStudent") + 
    theme_light() +
    theme(axis.text.x = element_text(angle =30, hjust = 1,size = 10), 
          axis.text.y = element_text(colour="grey",angle = 90, size = 12),
          axis.title.y = element_text(angle = 90, size = 15)) +
    theme(legend.position = "right")+
    theme(axis.title=element_text(size=20),
          axis.text = element_text(face = "bold",size = 16),
          axis.ticks.length=unit(.4, "cm"),
          axis.ticks = element_line(colour = "black", size = 1),
          #panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1.5),
          plot.margin = margin(1, 1, 1, 1, "cm"))
  
  c+stat_compare_means(method = 'wilcox.test')
  ggsave(paste0("/home/data/gaoyuzhen/Projects/ImmunePro/Figures_PHD/","full_TMB_LM.pdf"),height=6.5,width=6.5)
  
 # c+stat_compare_means(method = 't.test', 
                       #label.y = c(2.5,3,3.5),
                       #label = "p.signif",
                       #symnum.args=list(),
                      # comparisons = list(c("liver", "bone_marrow"),
                                       #   c("liver", "brain"),
                                        #  c("liver", "lung"),
                                         # c("liver", "skin")))
  #ggsave(file.path(LiverFeatures.path ,paste0("Figure1/",gene,".pdf")),height=6.5,width=6.5)
  ggsave(paste0("Figure1/TMB/","full_TMB_LM.pdf"),height=6.5,width=6.5)
  ggsave(paste0("Figure1/TMB/","full_TMB_LM.png"),height=6.5,width=6.5)
  

  
  ##msi
  colnames(tmp)
  tmp[,"fga"]<-as.numeric(tmp[,"fga"])
  topvealues<-median(tmp[,"fga"])+5*sd(tmp[,"fga"])
  tmp[,"fga"]<-ifelse(tmp[,"fga"]>  topvealues,  topvealues,tmp[,"fga"])
  tmp$MetaLiver<-as.factor(tmp$MetaLiver)
  ################
  c<- ggplot(tmp,
             aes_string(x="MetaLiver", y="fga", 
                        fill ="MetaLiver", 
                        color ="MetaLiver")) +
    geom_point(aes(shape=Cancer_Type), 
               position = position_jitter(width = .15,height=-0.7),
               size=1)+
    xlab("")+
    scale_shape_manual(values=shapes) +
    scale_color_manual(values=c("black","#E31A1C")) +
    scale_fill_manual(values=c("black","#E31A1C")) +
    geom_boxplot(notch = F, alpha = 0.65, width=0.2, outlier.size = 0.65)+
    geom_violin(position = position_nudge(x = 0),width=0.9,alpha = 0.4)+
    #ggtitle("CPSS score", "stratified by InternStudent") + 
    theme_light() +
    theme(axis.text.x = element_text(angle =30, hjust = 1,size = 10), 
          axis.text.y = element_text(colour="grey",angle = 90, size = 12),
          axis.title.y = element_text(angle = 90, size = 15)) +
    theme(legend.position = "right")+
    theme(axis.title=element_text(size=20),
          axis.text = element_text(face = "bold",size = 16),
          axis.ticks.length=unit(.4, "cm"),
          axis.ticks = element_line(colour = "black", size = 1),
          #panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1.5),
          plot.margin = margin(1, 1, 1, 1, "cm"))
  
  c+stat_compare_means(method = 'wilcox.test')
  
  ggsave(paste0("/home/data/gaoyuzhen/Projects/ImmunePro/Figures_PHD/","full_fga_LM.pdf"),height=6.5,width=6.5)
  
  
  ggsave(paste0("Figure1/TMB/","fga_LM.pdf"),height=6.5,width=6.5)
  ggsave(paste0("Figure1/TMB/","fga_LM.png"),height=6.5,width=6.5)
  
#################################################################################  
  ##msi
  colnames(tmp)
  tmp[,"msi_score"]<-as.numeric(tmp[,"msi_score"])
  topvealues<-median(tmp[,"msi_score"])+2*sd(tmp[,"msi_score"])
  tmp[,"msi_score"]<-ifelse(tmp[,"msi_score"]>  topvealues,  topvealues,tmp[,"msi_score"])
  tmp$MetaLiver<-as.factor(tmp$MetaLiver)
  ################
  c<- ggplot(tmp,
             aes_string(x="MetaLiver", y="msi_score", 
                        fill ="MetaLiver", 
                        color ="MetaLiver")) +
    geom_point(aes(shape=Cancer_Type), 
               position = position_jitter(width = .15,height=-0.7),
               size=1)+
    xlab("")+
    scale_shape_manual(values=shapes) +
    scale_color_manual(values=c(mycol)) +
    scale_fill_manual(values=c(mycol)) +
    geom_boxplot(notch = F, alpha = 0.65, width=0.2, outlier.size = 0.65)+
    geom_violin(position = position_nudge(x = 0),width=0.9,alpha = 0.4)+
    #ggtitle("CPSS score", "stratified by InternStudent") + 
    theme_light() +
    theme(axis.text.x = element_text(angle =30, hjust = 1,size = 10), 
          axis.text.y = element_text(colour="grey",angle = 90, size = 12),
          axis.title.y = element_text(angle = 90, size = 15)) +
    theme(legend.position = "right")+
    theme(axis.title=element_text(size=20),
          axis.text = element_text(face = "bold",size = 16),
          axis.ticks.length=unit(.4, "cm"),
          axis.ticks = element_line(colour = "black", size = 1),
          #panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1.5),
          plot.margin = margin(1, 1, 1, 1, "cm"))
  
  c+stat_compare_means(method = 'wilcox.test')
  ggsave(paste0("Figure1/TMB/","msi_score_LM.pdf"),height=6.5,width=6.5)
  ggsave(paste0("Figure1/TMB/","msi_score_LM.png"),height=6.5,width=6.5)
  