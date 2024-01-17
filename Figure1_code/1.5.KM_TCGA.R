#####
library(survminer)
library(survival)
library(readxl)
library(ggplot2)
library(tibble)
library(scales)
library(ggrepel)
library(forcats)
library(reshape2)
library(cowplot)
library(data.table)
library(tidyr)
library(readr)

TCGA<-read_excel("Figure1/PANCANCER_information.xlsx")
colnames(TCGA)
TCGA$OS.time<-as.numeric(TCGA$OS.time)
TCGA$OS<-ifelse(TCGA$OS=="Dead",1,0)
TCGA$MetaLiver<-ifelse(TCGA$A8_New_Event_Tissue=="Liver",1,0)

mycol<-c('#D95F02','#7570B3','#E7298A', '#66A61E','#E6AB02','#A6761D','#666666','#377EB8',
         '#4DAF4A', '#984EA3','#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999','#8DD3C7',
         '#FFFFB3','#BEBADA', '#FB8072','#80B1D3','#FDB462','#B3DE69', '#FCCDE5',
         '#D9D9D9','#BC80BD','#CCEBC5', '#FFED6F','#66C2A5', '#FC8D62','#8DA0CB','#E78AC3',
         '#A6D854','#FFD92F', '#E5C494', '#B3B3B3','#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C',
         '#1B9E77','#7FC97F','#BEAED4','#FDC086','#FFFF99','#386CB0',
         '#FB9A99', '#E31A1C','#FDBF6F', '#FF7F00', '#CAB2D6', '#6A3D9A', '#FFFF99', '#B15928')
show_col(mycol)
#######################################
Clinicloutcome<-"Figure1"
data_unimmune<-TCGA
freq<-data.frame(table(data_unimmune$Type))
data_unimmune<-data_unimmune[data_unimmune$Type %in% freq$Var1[freq$Freq>11],]
#cancernames<-c("bladder cancer","breast cancer","colorectal cancer", 'esophagogastric cancer','head and neck cancer',"melanoma", "non-small cell lung cancer",'renal cell carcinoma')
cox.sig<-c()
set<-"all"
for (set in unique(data_unimmune$Type)){
  data_immunotherapy<-data_unimmune[data_unimmune$Type==set,]
  data_immunotherapy<-data_unimmune
  Gcox1<-coxph(Surv(OS.time,OS)~MetaLiver,data=data_immunotherapy)
  GSum<-summary(Gcox1)
  HR<-round(GSum$coefficients[,2],3)
  Pvalue<-round(GSum$coefficients[,5],3)
  CI<-paste0(round(GSum$conf.int[,3:4],3),collapse='-')
  coeff<-round(GSum$coefficients[1],3)
  se<-GSum$coefficients[1,3]
  low<-round(GSum$conf.int[,3],3)
  up<-round(GSum$conf.int[,4],3)
  cox.p<-data.frame(#'dataset'= set,
    'Hazard Ratio'= HR,
    'CI95'=CI,
    "coeff"=coeff,
    "se"=se,
    "low"=low,
    "up"=up,
    'P-value'=Pvalue,
    "set"=set,
    "nofliver"=data.frame(table(data_immunotherapy$MetaLiver))[2,2],
    "prec"=prop.table(table(data_immunotherapy$MetaLiver))[2],
    "numberofpatients"=nrow(data_immunotherapy))
  cox.p
  cox.sig=rbind(cox.sig,cox.p)
  message(set,"cox regression is Done!")
  ##plot the km curves
  kmfit<- survfit(Surv(OS.time/30,OS)~MetaLiver,data=data_immunotherapy)
  p.val<-round(Pvalue,3)
  HR <- paste("Hazard Ratio= ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low,2), round(up,2), sep = " - "), sep = "")
  #if (!dir.exists("metastasis_COHORT")){
  #  dir.create("./metastasis_COHORT")
  # }
  #if(p.val<0.05){
  pdf(file.path(Clinicloutcome,paste0(set,"_Unimmune_OS.pdf")),width = 6.49, height = 6.58,onefile = F )
  print(ggsurvplot(kmfit,#surv.median.line = "hv", # Add medians survival
                   # Change legends: title & labels
                   main = "Survival curve",
                   #legend.title = paste0(set,"-",sig),
                   legend.labs = c("NoLiverMeta","LiverMeta"),
                   #xlim = c(0, 2000),
                   xlab="Time(months)",
                   ylab="Overall Survival",
                   size = 1,
                   #fun="cumhaz",
                   #fun='event',
                   # Add p-value and tervals
                   pval = TRUE,
                   #test.for.trend = TRUE,###group more than 2groups
                   break.time.by = 30,
                   #conf.int = TRUE,
                   #group.by=,
                   # Add risk table
                   risk.table = TRUE,
                   tables.height = 0.185,
                   tables.theme = theme_cleantable(),
                   palette = c("black","#E31A1C"),
                   #ggtheme = theme_bw(), # Change ggplot2 theme
                   #font.title="OS",
                   font.main =15,
                   font.x =  15,
                   font.y = 15,
                   font.tickslab =15
                   #在左下???标出pvalue、HR???95% CI
                   #???小的p value标为p < 0.001
                   #pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                   #paste("p = ",round(p.val,3), sep = "")), HR, CI, sep = "\n")
  ) )
  invisible(dev.off())
}
###integrate the data of them.
write.csv(cox.sig,file="Figure1/cox.sig_TCGA.csv")

pie(table(data_unimmune$MetaLiver))

pie( table(rt$Type))
#### single cohort analysis for cox regression

