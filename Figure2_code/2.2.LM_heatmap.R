###################
###无监督聚类for
load("~/Projects/ImmuneProGPR/pathways_for_gsva_mets.Rdata")
load("~/Projects/MET500Projects/phe_met.Rdata")
load("~/Projects/MET500Projects/metExp500.Rdata")
load("~/Projects/MET500Projects/immunefeaturs.Rdata")
MET500_geneExpression_M_meta_plus <- read.delim2("~/Projects/MET500Projects/Clinicaldata/MET500_geneExpression_M.meta.plus.txt")
load("Data/data_LM.Rdata")
#################
#data_LM<-read.csv("Figure2.1/export_phe_figure2.1.csv",row.names = 1)
#data_LM_immunefeatures<-immunefeaturs[immunefeaturs$Sample_id %in% rownames(data_LM),]### selected patients for diff
#save(data_LM,data_LM_immunefeatures,file = "Data/data_LM.Rdata")

#data_LM_immunefeatures<-data_LM_immunefeatures[,-c(233:258)]
#phe_met$Sample_id<-gsub("-",".",phe_met$Sample_id)
#data_LM_immunefeatures<-merge(data_LM_immunefeatures,phe_met,by="Sample_id")


library(limma)
rt<-data_LM_immunefeatures
rt$score_group<-ifelse(rt$biopsy_tissue_bio=='liver',"yes","no")


group_list<-rt[,"score_group"]
design <- model.matrix(~0+factor(group_list))#�ѷ�����������ʽ
colnames(design)=levels(factor(group_list))#��design����Ϊ������Ϣ


exprSet<-t(rt[,c(2:93)])
colnames(exprSet)<-rt$Sample_id

exprSet<-data.frame(exprSet)
exprSet[]<-apply(exprSet,2,as.numeric)
rownames(design)=colnames(exprSet)#
contrast.matrix<-makeContrasts("yes-no",levels = design)
#contrast.matrix #
deg = function(exprSet,design,contrast.matrix){#����һ����deg�ĺ���
  ##step1#logFC����ǰ�Ķ��ٱ�,���в��컯�Ƚ�
  fit <- lmFit(exprSet,design)
  ##step2
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  ##��һ������Ҫ����ҿ������п���Ч��
  
  fit2 <- eBayes(fit2)  ## default no trend !!!
  ##eBayes() with trend=TRUE
  ##step3
  tempOutput = topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput) 
  #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
  head(nrDEG)
  return(nrDEG)
}
deg = deg(exprSet,design,contrast.matrix)
##
library(tidyverse)
library(reshape)
library(ggplot2)



celltype<-read.delim2("/home/data/gaoyuzhen/Projects/ImmuneProGPR/Celltype.txt")
write.csv(outTab,file = "outTab.csv")
outTab<-read.csv("outTab.csv",row.names = 1)
names(df)
table(outTab$type)
df<-outTab[outTab$type=="Lymphoid",]
df$genes<-rownames(df)
df$logp.value<--log(df$adj.P.Val+0.0001)
df$lm<-ifelse(df$logFC>0,"1","0")
p<-ggplot(df, aes(x=genes, y=logp.value, fill=lm,colour=lm)) +
   geom_bar(stat="identity", linetype=ifelse(df$adj.P.Val < 0.001, "solid", "dashed")) +
  scale_color_manual(values=c("black","#E31A1C")) +
  scale_fill_manual(values=c("black","#E31A1C")) +
  #geom_histogram(binwidth=1, stat='identity') +
  theme_light() + 
  
  #scale_fill_gradient(low='white', high='red', limits=c(0,10)) + 
  theme(axis.title.y=element_text(angle=0))

p + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+geom_hline(yintercept = -log10(pvalCut), color="grey40", 
                                                                               linetype="longdash", lwd = 0.5)
###################
table(outTab$type)
df<-outTab[outTab$type=="Lymphoid",]


df<-outTab[outTab$type!="Lymphoid" & outTab$type!="Score",]
df$genes<-rownames(df)
df$logp.value<--log(df$adj.P.Val+0.0001)
df$lm<-ifelse(df$logFC>0,"1","0")
p<-ggplot(df, aes(x=genes, y=logp.value, fill=lm,colour=lm)) +
  geom_bar(stat="identity", linetype=ifelse(df$adj.P.Val < 0.001, "solid", "dashed")) +
  scale_color_manual(values=c("black","#E31A1C")) +
  scale_fill_manual(values=c("black","#E31A1C")) +
  #geom_histogram(binwidth=1, stat='identity') +
  theme_light() + 
  
  #scale_fill_gradient(low='white', high='red', limits=c(0,10)) + 
  theme(axis.title.y=element_text(angle=0))

p + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))+geom_hline(yintercept = -log10(pvalCut), color="grey40", 
                                                                               linetype="longdash", lwd = 0.5)


########################################
p + coord_polar()
p + coord_polar() + aes(x=reorder(genes, logp.value)) +
  theme(axis.text.x = element_text(angle=-20)) 
ggsave("Lymphoid.pdf")
ggsave("NoLymphoid.pdf")

###################################################


##########
library(ggplot2)
df$logFC<-as.character(as.numeric(factor(df$logFC)))
df$logFC<-as.numeric(df$logFC)
df$P.Value<-as.numeric(df$P.Value)
######special volcano plot
rownames(df)
x<-df
x$label<- rownames(x)
head(x)

#plot_mode <- "classic" #?????
plot_mode <- "advanced" #?????
#
mid<-40
logFCcut <- 75 #log2-foldchange
pvalCut <- 0.05 #P.value
adjPcut <- 0.05 #adj.P.value

#for advanced mode
logFCcut2 <- 75
pvalCut2 <- 0.01
pvalCut3 <- 0.001

#########
xmin <- (range(x$logFC)[1]- (range(x$logFC)[1]+ 10))
xmax <- (range(x$logFC)[1]+ (10-range(x$logFC)[1]))
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
  cols[x$P.Value < pvalCut & x$logFC >mid]<- "#FB9A99"
  cols[x$P.Value < pvalCut2 & x$logFC > mid]<- "#ED4F4F"
  cols[x$P.Value < pvalCut & x$logFC < mid]<- "#B2DF8A"
  cols[x$P.Value < pvalCut2 & x$logFC < mid]<- "#329E3F"
  color_transparent <- adjustcolor(cols, alpha.f = 0.5)
  x$color_transparent <- color_transparent
  
  # ?}?s???setting for size
  n1 <- length(x[, 1])
  size <- rep(1, n1)
  
  #?????????????????
  size[x$P.Value < pvalCut & x$logFC > mid]<- 2
  size[x$P.Value < pvalCut2 & x$logFC > mid]<- 4
  size[x$P.Value < pvalCut3 & x$logFC >mid]<- 6
  size[x$P.Value < pvalCut & x$logFC < mid]<- 2
  size[x$P.Value < pvalCut2 & x$logFC < mid]<- 4
  size[x$P.Value < pvalCut3 & x$logFC <mid]<- 6
  
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
    geom_vline(xintercept = c(mid, mid), color="grey40", 
               linetype="longdash", lwd = 0.5) +
    geom_hline(yintercept = -log10(pvalCut), color="grey40", 
               linetype="longdash", lwd = 0.5)+
    geom_hline(yintercept = -log10(pvalCut2), color="grey40", 
               linetype="longdash", lwd = 0.5)+
    geom_hline(yintercept = -log10(pvalCut3), color="grey40", 
               linetype="longdash", lwd = 0.5)
}
p1

##########interation of genes in the BCC tumor and adjucent tissue
library(ggplot2)
library(ggrepel)
library(ggthemes)
n = 75
p1 + geom_text_repel(data=x,max.overlaps=6,aes(x = logFC, y = -log10(P.Value), 
                                label = ifelse(-log(x$P.Value)>3, rownames(x),"")),
                     colour="darkred", size = 5, box.padding = unit(0.35, "lines"), 
                     point.padding = unit(0.3, "lines"))

ggsave(file = "Lymphoid_LM_cells.pdf",width =6,height =6)####

ggsave(file = "Non_Lymphoid_LM_cells.pdf",width =6,height =6)####
