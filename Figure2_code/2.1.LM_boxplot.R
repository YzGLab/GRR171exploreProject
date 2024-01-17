####
load("~/Projects/ImmuneProGPR/pathways_for_gsva_mets.Rdata")
load("~/Projects/MET500Projects/phe_met.Rdata")
load("~/Projects/MET500Projects/metExp500.Rdata")
load("~/Projects/MET500Projects/immunefeaturs.Rdata")
#load("~/Projects/MET500Projects/immunefeaturs_GEO.Rdata")
###for cancer pathway
rt2<-immunefeaturs
phe_met<-phe_met[phe_met$tissue %in% c("breast","colon","lung","pancreas","prostate","esophagus","stomach","skin"),]
phe_met$Sample_id<-gsub("-",".",phe_met$Sample_id)
rownames( phe_met)<- phe_met$Sample_id
#phe_met<-phe_met[phe_met$tc>0.6,]
phe_met<-phe_met[rownames(phe_met) %in% rownames(rt2),]
rt2<-rt2[rownames(rt2) %in% phe_met$Sample_id,]
##################################################################
dim(phe_met)
############
type1<-phe_met[,c(13,12,7)]
rt2<-data.frame(rt2)
rt2<-cbind(rt2,type1)

shape_level <- nlevels(factor(rt2$cohort))
if (shape_level < 15){
  shapes = (0:shape_level) %% 15
} else{
  shapes = c(0:14,c((15:shape_level) %% 110 + 18))
}


comparVariables<-colnames(rt2)[c(193:207)]
for(gene in comparVariables){
  print(gene)
  rt2[,gene]<-as.numeric(rt2[,gene])
  topvealues<-median(rt2[,gene])+4*sd(rt2[,gene])
  rt2[,gene]<-ifelse(rt2[,gene]>  topvealues,  topvealues,rt2[,gene])
  ################
  c<- ggplot(rt2,
             aes_string(x="biopsy_tissue_bio", y=gene, 
                        fill ="biopsy_tissue_bio", 
                        color ="biopsy_tissue_bio")) +
    geom_jitter(aes_string(shape="cohort"), 
                position = position_jitter(width = .15,height=-0.7),
                size=1)+
    xlab("")+
    scale_shape_manual(values=shapes) +
    scale_color_manual(values=c("#E31A1C","black")) +
    scale_fill_manual(values=c("#E31A1C","black")) +
    geom_boxplot(notch = F, alpha = 0.65, width=0.2, outlier.size = 0.85)+
    geom_violin(position = position_nudge(x = 0),width=0.9,alpha = 0.4)+
    #scale_y_continuous(limits = quantile(rt[,gene], c(0.05, 0.95)))+
    #ggtitle("CPSS score", "stratified by InternStudent") + 
    theme_classic() +
    theme(axis.text.x = element_text(angle =30, hjust = 1,size = 10), 
          axis.text.y = element_text(colour="grey",angle = 90, size = 12),
          axis.title.y = element_text(angle = 90, size = 15)) +
    #theme(legend.position = "top")
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
  c+stat_compare_means(method = 't.test')
  ggsave(file.path(LiverFeatures.path,paste0("Pan_LM_",gene,".pdf")),height=5.5,width=6.24)
}
######

comparVariables<-colnames(rt2)[c(258:369)]
for(gene in comparVariables){
  print(gene)
  rt2[,gene]<-as.numeric(rt2[,gene])
  topvealues<-median(rt2[,gene])+4*sd(rt2[,gene])
  rt2[,gene]<-ifelse(rt2[,gene]>  topvealues,  topvealues,rt2[,gene])
  ################
  c<- ggplot(rt2,
             aes_string(x="biopsy_tissue_bio", y=gene, 
                        fill ="biopsy_tissue_bio", 
                        color ="biopsy_tissue_bio")) +
    geom_jitter(aes_string(shape="cohort"), 
                position = position_jitter(width = .15,height=-0.7),
                size=1)+
    xlab("")+
    scale_shape_manual(values=shapes) +
    scale_color_manual(values=c("#E31A1C","black")) +
    scale_fill_manual(values=c("#E31A1C","black")) +
    geom_boxplot(notch = F, alpha = 0.65, width=0.2, outlier.size = 0.85)+
    geom_violin(position = position_nudge(x = 0),width=0.9,alpha = 0.4)+
    #scale_y_continuous(limits = quantile(rt[,gene], c(0.05, 0.95)))+
    #ggtitle("CPSS score", "stratified by InternStudent") + 
    theme_classic() +
    theme(axis.text.x = element_text(angle =30, hjust = 1,size = 10), 
          axis.text.y = element_text(colour="grey",angle = 90, size = 12),
          axis.title.y = element_text(angle = 90, size = 15)) +
    #theme(legend.position = "top")
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
  c+stat_compare_means(method = 't.test')
  ggsave(paste0("Figure2/metabolism/","Pan_LM_",gene,".pdf"),height=5.5,width=6.24)
}
##
###
comparVariables<-colnames(rt2)[c(142)]
for(gene in comparVariables){
  print(gene)
  rt2[,gene]<-as.numeric(rt2[,gene])
  topvealues<-median(rt2[,gene])+3*sd(rt2[,gene])
  rt2[,gene]<-ifelse(rt2[,gene]>  topvealues,  topvealues,rt2[,gene])
  ################
  c<- ggplot(rt2,
             aes_string(x="biopsy_tissue_bio", y=gene, 
                        fill ="biopsy_tissue_bio", 
                        color ="biopsy_tissue_bio")) +
    geom_jitter(aes_string(shape="cohort"), 
                position = position_jitter(width = .15,height=-0.7),
                size=1)+
    xlab("")+
    scale_shape_manual(values=shapes) +
    scale_color_manual(values=c("#E31A1C","black")) +
    scale_fill_manual(values=c("#E31A1C","black")) +
    geom_boxplot(notch = F, alpha = 0.65, width=0.2, outlier.size = 0.85)+
    geom_violin(position = position_nudge(x = 0),width=0.9,alpha = 0.4)+
    #scale_y_continuous(limits = quantile(rt[,gene], c(0.05, 0.95)))+
    #ggtitle("CPSS score", "stratified by InternStudent") + 
    theme_classic() +
    theme(axis.text.x = element_text(angle =30, hjust = 1,size = 10), 
          axis.text.y = element_text(colour="grey",angle = 90, size = 12),
          axis.title.y = element_text(angle = 90, size = 15)) +
    #theme(legend.position = "top")
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
  c+stat_compare_means(method = 't.test')
  ggsave(paste0("Figure2/Immunecellsxcell/","Pan_LM_",gene,".pdf"),height=5.5,width=6.24)
  ggsave(paste0("Figure2/Immunecellsxcell/","Pan_LM_",gene,".png"),height=5.5,width=6.24)
}


###
comparVariables<-colnames(rt2)[c(2:11)]
for(gene in comparVariables){
  print(gene)
  rt2[,gene]<-as.numeric(rt2[,gene])
  topvealues<-median(rt2[,gene])+3*sd(rt2[,gene])
  rt2[,gene]<-ifelse(rt2[,gene]>  topvealues,  topvealues,rt2[,gene])
  ################
  c<- ggplot(rt2,
             aes_string(x="biopsy_tissue_bio", y=gene, 
                        fill ="biopsy_tissue_bio", 
                        color ="biopsy_tissue_bio")) +
    geom_jitter(aes_string(shape="cohort"), 
                position = position_jitter(width = .15,height=-0.7),
                size=1)+
    xlab("")+
    scale_shape_manual(values=shapes) +
    scale_color_manual(values=c("#E31A1C","black")) +
    scale_fill_manual(values=c("#E31A1C","black")) +
    geom_boxplot(notch = F, alpha = 0.65, width=0.2, outlier.size = 0.85)+
    geom_violin(position = position_nudge(x = 0),width=0.9,alpha = 0.4)+
    #scale_y_continuous(limits = quantile(rt[,gene], c(0.05, 0.95)))+
    #ggtitle("CPSS score", "stratified by InternStudent") + 
    theme_classic() +
    theme(axis.text.x = element_text(angle =30, hjust = 1,size = 10), 
          axis.text.y = element_text(colour="grey",angle = 90, size = 12),
          axis.title.y = element_text(angle = 90, size = 15)) +
    #theme(legend.position = "top")
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
  c+stat_compare_means(method = 't.test')
  ggsave(paste0("Figure2/Immunecells/","Pan_LM_",gene,".pdf"),height=5.5,width=6.24)
}

###
comparVariables<-colnames(rt2)[c(13:78)]
for(gene in comparVariables){
  print(gene)
  rt2[,gene]<-as.numeric(rt2[,gene])
  topvealues<-median(rt2[,gene])+3*sd(rt2[,gene])
  rt2[,gene]<-ifelse(rt2[,gene]>  topvealues,  topvealues,rt2[,gene])
  ################
  c<- ggplot(rt2,
             aes_string(x="biopsy_tissue_bio", y=gene, 
                        fill ="biopsy_tissue_bio", 
                        color ="biopsy_tissue_bio")) +
    geom_jitter(aes_string(shape="cohort"), 
                position = position_jitter(width = .15,height=-0.7),
                size=1)+
    xlab("")+
    scale_shape_manual(values=shapes) +
    scale_color_manual(values=c("#E31A1C","black")) +
    scale_fill_manual(values=c("#E31A1C","black")) +
    geom_boxplot(notch = F, alpha = 0.65, width=0.2, outlier.size = 0.85)+
    geom_violin(position = position_nudge(x = 0),width=0.9,alpha = 0.4)+
    #scale_y_continuous(limits = quantile(rt[,gene], c(0.05, 0.95)))+
    #ggtitle("CPSS score", "stratified by InternStudent") + 
    theme_classic() +
    theme(axis.text.x = element_text(angle =30, hjust = 1,size = 10), 
          axis.text.y = element_text(colour="grey",angle = 90, size = 12),
          axis.title.y = element_text(angle = 90, size = 15)) +
    #theme(legend.position = "top")
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
  c+stat_compare_means(method = 't.test')
  ggsave(paste0("Figure2/Immunecellsxcell/","Pan_LM_",gene,".pdf"),height=5.5,width=6.24)
  ggsave(paste0("Figure2/Immunecellsxcell/","Pan_LM_",gene,".png"),height=5.5,width=6.24)
}


