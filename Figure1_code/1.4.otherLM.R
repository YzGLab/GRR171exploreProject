#####
##
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
mycol<-c('#D95F02','#7570B3','#E7298A', '#66A61E','#E6AB02','#A6761D','#666666','#377EB8',
         '#4DAF4A', '#984EA3','#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999','#8DD3C7',
         '#FFFFB3','#BEBADA', '#FB8072','#80B1D3','#FDB462','#B3DE69', '#FCCDE5',
         '#D9D9D9','#BC80BD','#CCEBC5', '#FFED6F','#66C2A5', '#FC8D62','#8DA0CB','#E78AC3',
         '#A6D854','#FFD92F', '#E5C494', '#B3B3B3','#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C',
         '#1B9E77','#7FC97F','#BEAED4','#FDC086','#FFFF99','#386CB0',
         '#FB9A99', '#E31A1C','#FDBF6F', '#FF7F00', '#CAB2D6', '#6A3D9A', '#FFFF99', '#B15928')
show_col(mycol)
#################
load("livemeta_cohort.Rdata")

Clinicloutcome<-"Figure1"

set<-names(livemeta_phe)[2]
sub_rt<-livemeta_phe[[set]]
sub_rt$OS.event<-sub_rt$OS.Event
livemeta_phe[[set]]<-sub_rt



library(readxl)
CRC_immun<- read_excel("Figure1/In_house_immunotheapy.xlsx",1)
CRC_immun<-CRC_immun[-c(1:9),]
CRC_no_immun<- read_excel("Figure1/In_house_immunotheapy.xlsx",2)
lung_immun<- read_excel("Figure1/In_house_immunotheapy.xlsx",3)
livemeta_phe[["inhouseCRC_immuno"]]<-CRC_immun
livemeta_phe[["inhouseCRC_non_immuno"]]<-CRC_no_immun
livemeta_phe[["inhouseLung_immuno"]]<-lung_immun
save(livemeta_phe,livemeta_cohort,file = "livemeta_cohort.Rdata")
#######################################

library(readxl)
library(readr)
data_immunotherapy_total<-read_excel("/home/data/gaoyuzhen/Projects/ImmunePro/Figure1/In_house_immunotheapy.xlsx",5)
table(data_immunotherapy_total$cancer_type)
#cancernames<-c("bladder cancer","breast cancer","colorectal cancer", 'esophagogastric cancer','head and neck cancer',"melanoma", "non-small cell lung cancer",'renal cell carcinoma')
cox.sig<-c()
for (set in unique(data_immunotherapy_total$cancer_type)){
data_immunotherapy<-data_immunotherapy_total[data_immunotherapy_total$cancer_type==set,]
Gcox1<-coxph(Surv(os_days,os_status)~MetaLiver,data=data_immunotherapy)
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
kmfit<- survfit(Surv(os_days,os_status)~MetaLiver,data=data_immunotherapy)
p.val<-round(Pvalue,3)
HR <- paste("Hazard Ratio= ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low,2), round(up,2), sep = " - "), sep = "")
#if (!dir.exists("metastasis_COHORT")){
#  dir.create("./metastasis_COHORT")
# }
#if(p.val<0.05){
pdf(file.path(Clinicloutcome,paste0(set,"-OS.pdf")),width = 6.49, height = 6.58,onefile = F )
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
                 break.time.by = 10,
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
data_immunotherapy_total<-data_immunotherapy_total[data_immunotherapy_total$cancer_type!="colorectal cancer_Sir_noimm",]
data_immunotherapy_total$cancer_type<-ifelse(data_immunotherapy_total$cancer_type %like% "colo","colorectal cancer",data_immunotherapy_total$cancer_type)
#data_immunotherapy_total$cancer_type<-ifelse(data_immunotherapy_total$cancer_type %like% "bladder","bladder cancer",data_immunotherapy_total$cancer_type)
data_immunotherapy_total$cancer_type<-ifelse(data_immunotherapy_total$cancer_type %like% "elanoma","melanoma",data_immunotherapy_total$cancer_type)
data_immunotherapy_total$cancer_type<-ifelse(data_immunotherapy_total$cancer_type %like% "non-small cell lung cancer","non-small cell lung cancer",data_immunotherapy_total$cancer_type)
data_immunotherapy_total$cancer_type<-ifelse(data_immunotherapy_total$cancer_type %like% "Braun_RCC","renal cell carcinoma",data_immunotherapy_total$cancer_type)
####################################

pie( table(data_immunotherapy_total$MetaLiver,data_immunotherapy_total$cancer_type))
pie( table(rt$cancer_type))
#### single cohort analysis for cox regression

