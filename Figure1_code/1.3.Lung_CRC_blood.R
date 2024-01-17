#
library(plyr)
library(survival)
library(readxl)

tmp<- read_excel("~/Projects/ImmunePro/Figure1/In_house_immunotheapy.xlsx",1) %>% data.frame()

#
colnames(tmp)

tmp[,"后CRP"]<-as.numeric(tmp[,"后CRP"])
#topvealues<-median(tmp[,"msi_score"])+5*sd(tmp[,"msi_score"])
#tmp[,"msi_score"]<-ifelse(tmp[,"msi_score"]>  topvealues,  topvealues,tmp[,"msi_score"])
tmp$MetaLiver<-ifelse(tmp$Livermeta=="Yes",1,0)
tmp$MetaLiver<-as.factor(tmp$MetaLiver)
################
c<- ggplot(tmp,
           aes_string(x="MetaLiver", y="后CEA", 
                      fill ="MetaLiver", 
                      color ="MetaLiver")) +
  geom_point(#aes(shape=Cancer_Type), 
             position = position_jitter(width = .15,height=-0.7),
             size=1)+
  xlab("")+
  ylab('CEA after immunotherapy')+
  #scale_shape_manual(values=shapes) +
  scale_color_manual(values= c("black","#E31A1C")) +
  scale_fill_manual(values= c("black","#E31A1C")) +
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
ggsave(paste0("/home/data/gaoyuzhen/Projects/ImmunePro/Figures_PHD/","CRC_后CEA_LM.pdf"),height=6.5,width=6.5)
##

################
c<- ggplot(tmp,
           aes_string(x="MetaLiver", y="后CRP", 
                      fill ="MetaLiver", 
                      color ="MetaLiver")) +
  geom_point(#aes(shape=Cancer_Type), 
    position = position_jitter(width = .15,height=-0.7),
    size=1)+
  xlab("")+
  ylab('CRP after immunotherapy')+
  #scale_shape_manual(values=shapes) +
  scale_color_manual(values= c("black","#E31A1C")) +
  scale_fill_manual(values= c("black","#E31A1C")) +
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
ggsave(paste0("/home/data/gaoyuzhen/Projects/ImmunePro/Figures_PHD/","CRC_后CRP_LM.pdf"),height=6.5,width=6.5)



tmp<- read_excel("~/Projects/ImmunePro/Figure1/In_house_immunotheapy.xlsx",3) %>% data.frame()

##

#
colnames(tmp)

tmp[,"LDH"]<-as.numeric(tmp[,"ldh"])
#topvealues<-median(tmp[,"msi_score"])+5*sd(tmp[,"msi_score"])
#tmp[,"msi_score"]<-ifelse(tmp[,"msi_score"]>  topvealues,  topvealues,tmp[,"msi_score"])
tmp$MetaLiver<-ifelse(tmp$Livermeta=="Yes",1,0)
tmp$MetaLiver<-as.factor(tmp$MetaLiver)
################
c<- ggplot(tmp,
           aes_string(x="MetaLiver", y="LDH", 
                      fill ="MetaLiver", 
                      color ="MetaLiver")) +
  geom_point(#aes(shape=Cancer_Type), 
    position = position_jitter(width = .15,height=-0.7),
    size=1)+
  xlab("")+
  ylab('LDH before immunotherapy')+
  #scale_shape_manual(values=shapes) +
  scale_color_manual(values= c("black","#E31A1C")) +
  scale_fill_manual(values= c("black","#E31A1C")) +
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
ggsave(paste0("/home/data/gaoyuzhen/Projects/ImmunePro/Figures_PHD/","Lung_LDH_LM.pdf"),height=6.5,width=6.5)
##

c<- ggplot(tmp,
           aes_string(x="MetaLiver", y="细胞角蛋白19片段", 
                      fill ="MetaLiver", 
                      color ="MetaLiver")) +
  geom_point(#aes(shape=Cancer_Type), 
    position = position_jitter(width = .15,height=-0.7),
    size=1)+
  xlab("")+
  ylab('CK19 before immunotherapy')+
  #scale_shape_manual(values=shapes) +
  scale_color_manual(values= c("black","#E31A1C")) +
  scale_fill_manual(values= c("black","#E31A1C")) +
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
ggsave(paste0("/home/data/gaoyuzhen/Projects/ImmunePro/Figures_PHD/","Lung_CK19_LM.pdf"),height=6.5,width=6.5)
