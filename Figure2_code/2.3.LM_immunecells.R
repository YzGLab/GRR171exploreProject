###establish the heatmap
library(circlize)
library(gghalves)
library(ggpubr)
library(RColorBrewer)
library(ComplexHeatmap)
library(readr)

load("~/Projects/MET500Projects/metExp500.Rdata")
load("~/Projects/MET500Projects/immunefeaturs.Rdata")
#load("~/Projects/MET500Projects/immunefeaturs_GEO.Rdata")
load("~/Projects/MET500Projects/phe_met.Rdata")
#############################################

rt2<-immunefeaturs
phe_met<-phe_met[phe_met$tissue %in% c("breast","colon","lung","pancreas","prostate","esophagus","stomach","skin"),]
phe_met$Sample_id<-gsub("-",".",phe_met$Sample_id)
rownames( phe_met)<- phe_met$Sample_id
#phe_met<-phe_met[phe_met$tc>0.6,]
phe_met<-phe_met[rownames(phe_met) %in% rownames(rt2),]
rt2<-rt2[rownames(rt2) %in% phe_met$Sample_id,]
rt2<-data.frame(rt2)
rt2<-cbind(rt2,type1)
#
shape_level <- nlevels(factor(rt2$cohort))
if (shape_level < 15){
  shapes = (0:shape_level) %% 15
} else{
  shapes = c(0:14,c((15:shape_level) %% 110 + 18))
}

comparVariables<-colnames(rt2)[c(13:78)]
p_s<-c()
##################################################
for(gene in comparVariables){
  print(gene)
  rt2[,gene]<-as.numeric(rt2[,gene])
  topvealues<-median(rt2[,gene])+3*sd(rt2[,gene])
  rt2[,gene]<-ifelse(rt2[,gene]>  topvealues,  topvealues,rt2[,gene])
  test<-t.test(rt2[,gene]~rt2$biopsy_tissue_bio)
  p<-c(mean(rt2[rt2$biopsy_tissue_bio=="liver",gene]),mean(rt2[rt2$biopsy_tissue_bio!="liver",gene]),test$p.value)
  p_s<-rbind(p_s,p)
  ################
  #ggsave(paste0("Figure2/Immunecellsxcell/","Pan_LM_",gene,".pdf"),height=5.5,width=6.24)
  #ggsave(paste0("Figure2/Immunecellsxcell/","Pan_LM_",gene,".png"),height=5.5,width=6.24)
}
####
rownames(p_s)<-comparVariables
p_s<-data.frame(p_s)
p_s[,c(1:3)]<-round(p_s[,c(1:3)],6)
p_s$level<-ifelse(p_s[,1]>p_s[,2],"LM","NonLM")
p_s[,1]<-as.numeric(p_s[,1])
p_s[,2]<-as.numeric(p_s[,2])
p_s$logFC<-log(p_s$X1/p_s$X2)
p_s$adj.P.value<-p.adjust(p_s$X3,method = "fdr")

sig.names<-rownames(p_s)[rownames(p_s) %in% rownames(p_s[p_s$V3<0.05,])]
p_s$pstar <- ifelse(p_s$adj.P.value < 0.05,
                    ifelse(p_s$adj.P.value < 0.01,paste0(rownames(p_s),"**"),paste0(rownames(p_s),"*")),
                    "a")
rownames(p_s)<-ifelse(p_s$pstar=="a",rownames(p_s),p_s$pstar)
outTab<-p_s
library(ggplot2)
colnames(outTab)
#rownames(outTab)<-outTab$gene
outTab$logFC<-as.numeric(outTab$logFC)
outTab$P.Value<-as.numeric(outTab$adj.P.value)
######special volcano plot
x<-outTab
x$label<- rownames(x)
head(x)

#plot_mode <- "classic" #?????
plot_mode <- "advanced" #?????
#
logFCcut <- 0.5 #log2-foldchange
pvalCut <- 0.05 #P.value
adjPcut <- 0.05 #adj.P.value

#for advanced mode
logFCcut2 <- 1
logFCcut3 <- 2
pvalCut2 <- 0.05
pvalCut3 <- 0.001

#??x??y?S?????????????????????
xmin <- (range(x$logFC)[1]- (range(x$logFC)[1]+ 10))
xmax <- (range(x$logFC)[1]+ (10-range(x$logFC)[1]))

#xmax <- 15
ymin <- 0
ymax <- -log10(x$P.Value)[3] * 1.1

mycol <- c("darkgreen","chocolate4","blueviolet","#223D6C","#D20A13",
           "#088247","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D",
           "#6E568C","#E0367A","#D8D155","#64495D","#7CC767")
if (plot_mode == "classic"){
  # ??????????setting for color
  x$color_transparent <- ifelse((x$P.Value < pvalCut & x$logFC > logFCcut), "red", ifelse((x$P.Value < pvalCut & x$logFC < -logFCcut), "blue","grey30"))
  # ??????????setting for size
  size <- ifelse((x$P.Value < pvalCut & abs(x$logFC) > logFCcut), 4, 2)
  
} else if (plot_mode == "advanced") {
  # ?}?s???setting for color
  n1 <- length(x[, 1])
  cols <- rep("grey30", n1)
  names(cols)<- rownames(x)
  
  #?????????????
  cols[x$P.Value < pvalCut & x$logFC >logFCcut]<- "#FB9A99"
  cols[x$P.Value < pvalCut2 & x$logFC > logFCcut2]<- "#ED4F4F"
  cols[x$P.Value < pvalCut & x$logFC < -logFCcut]<- "#B2DF8A"
  cols[x$P.Value < pvalCut2 & x$logFC < -logFCcut2]<- "#329E3F"
  color_transparent <- adjustcolor(cols, alpha.f = 0.5)
  x$color_transparent <- color_transparent
  
  # ?}?s???setting for size
  n1 <- length(x[, 1])
  size <- rep(1, n1)
  
  #?????????????????
  size[x$P.Value < pvalCut & x$logFC > logFCcut]<- 2
  size[x$P.Value < pvalCut2 & x$logFC > logFCcut2]<- 4
  size[x$P.Value < pvalCut3 & x$logFC > logFCcut3]<- 6
  size[x$P.Value < pvalCut & x$logFC < -logFCcut]<- 2
  size[x$P.Value < pvalCut2 & x$logFC < -logFCcut2]<- 4
  size[x$P.Value < pvalCut3 & x$logFC < -logFCcut3]<- 6
  
} else {
  stop("Unsupport mode")
}

# Construct the plot object
p1 <- ggplot(data=x, aes(logFC, -log10(P.Value), label = label)) +
  geom_point(alpha = 0.8, size = size, colour = x$color_transparent) +
  #labs(x=bquote(~Log[2]~"(fold change)"), y=bquote(~-Log[10]~italic("P-value")), title="") + 
  #ylim(c(ymin,ymax)) + 
  #xlim(c(-4,4))+
  #scale_x_continuous(
  # breaks = c(-5, -3, -logFCcut, 0, logFCcut, 3, 5), #??????????????
  #labels = c(-5, -3, -logFCcut, 0, logFCcut, 3, 5),
  #limits = c(-11, 11) #x???????????????????
  #) +
  #??????????????????
  #xlim(c(xmin, xmax)) + 
  
#??????????
#geom_vline(xintercept = c(-logFCcut, logFCcut), color="grey40", 
#    linetype="longdash", lwd = 0.5) + #????????????
# geom_hline(yintercept = -log10(pvalCut), color="grey40", 
#linetype="longdash", lwd = 0.5) +

theme_bw(base_size = 12#, base_family = "Times" #???????
) +
  theme(panel.grid=element_blank())
p1

if (plot_mode == "advanced") {
  p1 <- p1 + 
    geom_vline(xintercept = c(-logFCcut2, logFCcut2), color="grey40", 
               linetype="longdash", lwd = 0.5) +
    geom_hline(yintercept = -log10(pvalCut2), color="grey40", 
               linetype="longdash", lwd = 0.5)
}
p1

##########interation of genes in the BCC tumor and adjucent tissue
library(ggplot2)
library(ggrepel)
library(ggthemes)
n = 1
p1 + geom_text_repel(data=x,aes(x = logFC, y = -log10(P.Value), 
                                label = ifelse(abs(logFC) > n, rownames(x),"")),
                     colour="darkred", size = 5, box.padding = unit(0.35, "lines"), 
                     point.padding = unit(0.3, "lines"))

ggsave("diff_immunecells.pdf",width=9,height=9)
selectgenes<-data.frame(rownames(data))

###
targetgenes<<-x[x$logFC>1.5,]
#targetgenes$label
#targetgenes<-targetgenes[sort(targetgenes$logFC,decreasing = F),]


p2 <- p1 + 
  # ??????????????????????????
  geom_point(data = selectgenes, alpha = 1, size = 4.6, shape = 1, 
             stroke = 1, #????
             color = "black") +
  
  # ???????????????????
  scale_color_manual(values = mycol) + 
  geom_text_repel(data = selectgenes, 
                  show.legend = FALSE, #????????
                  size = 5, box.padding = unit(0.35, "lines"), 
                  point.padding = unit(0.3, "lines")) +
  guides(color=guide_legend(title = NULL)) 

p2




library(pheatmap)
tmp<-data.frame(p_s[,c(1,2)])
tmp<-t(scale(t(tmp)))
tmp<-tmp[tmp[,1]!="NaN",]
pdf("Figure2/immunecellxcellExpression.pdf",width = 6.15,height=11.69)
Heatmap(tmp,name = "Raltive AveExpression",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 12)
        # row_split = as$splitgroups ,
        #row_names_side =  "left",
        #right_annotation = rowAnn
)
dev.off()
####


###
comparVariables<-colnames(rt2)[c(2:11)]
comparVariables
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
    facet_wrap(~cohort)+
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
  ggsave(paste0("Figure2/Immunecells/","Single_Pan_LM_",gene,".pdf"),height=5.5,width=6.24)
}