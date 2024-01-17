###
library(plyr)
library(survival)
library(readxl)
CRC_immun<- read_excel("~/Projects/ImmunePro/Figure1/In_house_immunotheapy.xlsx",1)
################


citPie <- prop.table(table(CRC_immun$Livermeta))
print(citPie)
pie(citPie,
    col=c("black","#E31A1C"))


citPie <- prop.table(table(CRC_immun$LungMeta,CRC_immun$Livermeta))
print(citPie)
pdf(file.path(Clinicloutcome,"CRC_circle_OS.pdf"),width = 6, height = 6,onefile = F )
pie(citPie,border = TRUE,
    col=c("black","grey","#E31A1C","darkblue"))
invisible(dev.off())
####
citPie <- prop.table(table(CRC_immun$Livermeta))
print(citPie)
pdf(file.path(Clinicloutcome,"CRC_circle_liver_OS.pdf"),width = 6, height = 6,onefile = F )
pie(citPie,border = TRUE,
    col=c("black","grey","#E31A1C","darkblue"))
invisible(dev.off())




colnames(CRC_immun)
Gcox1<-coxph(Surv(OS,OS.event)~Livermeta+LungMeta,data=CRC_immun)
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
  'P-value'=Pvalue
)
cox.p

##plot the km curves
kmfit<- survfit(Surv(OS,OS.event)~Livermeta+LungMeta,data=CRC_immun)
p.val<-round(Pvalue,3)
HR <- paste("Hazard Ratio= ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low,2), round(up,2), sep = " - "), sep = "")
#if (!dir.exists("metastasis_COHORT")){
#  dir.create("./metastasis_COHORT")
# }
#if(p.val<0.05){
Clinicloutcome<-"Figure1"

pdf(file.path(Clinicloutcome,"CRC_Sir_OS.pdf"),width = 6.49, height = 6.58,onefile = F )
print(ggsurvplot(kmfit,#surv.median.line = "hv", # Add medians survival
                 # Change legends: title & labels
                 main = "Survival curve",
                 #legend.title = paste0(set,"-",sig),
                 legend.labs = c("NoLiver,NoLung","NoLiver,Lung","Liver,NoLung","Liver,Lung"),
                 #xlim = c(0, 2000),
                 xlab="Time(months)",
                 ylab="OS(%)",
                 size = 1,
                 #fun="cumhaz",
                 #fun='event',
                 # Add p-value and tervals
                 pval = TRUE,
                 #test.for.trend = TRUE,###group more than 2groups
                 break.time.by = 6,
                 #conf.int = TRUE,
                 #group.by=,
                 # Add risk table
                 risk.table = TRUE,
                 tables.height = 0.185,
                 tables.theme = theme_cleantable(),
                 palette = c("black","grey","#E31A1C","darkblue"),
                 #ggtheme = theme_bw(), # Change ggplot2 theme
                 #font.title="OS",
                 font.main =15,
                 font.x =  15,
                 font.y = 15,
                 font.tickslab =15,
                 #在左下???标出pvalue、HR???95% CI
                 #???小的p value标为p < 0.001
                 #pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                 # paste("p = ",round(p.val,3), sep = "")), HR, CI, sep = "\n")
) 
)
invisible(dev.off())
##
Gcox1<-coxph(Surv(OS,OS.event)~Livermeta,data=CRC_immun)
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
  'P-value'=Pvalue
)
cox.p

##plot the km curves
kmfit<- survfit(Surv(OS,OS.event)~Livermeta,data=CRC_immun)
p.val<-round(Pvalue,3)
HR <- paste("Hazard Ratio= ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low,2), round(up,2), sep = " - "), sep = "")
#if (!dir.exists("metastasis_COHORT")){
#  dir.create("./metastasis_COHORT")
# }
#if(p.val<0.05){
Clinicloutcome<-"Figure1"

pdf(file.path(Clinicloutcome,"CRC_Sir_OS——LM.pdf"),width = 6.49, height = 6.58,onefile = F )
print(ggsurvplot(kmfit,#surv.median.line = "hv", # Add medians survival
                 # Change legends: title & labels
                 main = "Survival curve",
                 #legend.title = paste0(set,"-",sig),
                 legend.labs = c("NoLiverMeta","LiverMeta"),
                 #xlim = c(0, 2000),
                 xlab="Time(months)",
                 ylab="OS(%)",
                 size = 1,
                 #fun="cumhaz",
                 #fun='event',
                 # Add p-value and tervals
                 pval = TRUE,
                 #test.for.trend = TRUE,###group more than 2groups
                 break.time.by = 6,
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
                 font.tickslab =15,
                 #在左下???标出pvalue、HR???95% CI
                 #???小的p value标为p < 0.001
                 #pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                 # paste("p = ",round(p.val,3), sep = "")), HR, CI, sep = "\n")
) 
)
invisible(dev.off())
###
Gcox1<-coxph(Surv(OS,OS.event)~LungMeta,data=CRC_immun)
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
  'P-value'=Pvalue
)
cox.p
kmfit<- survfit(Surv(OS,OS.event)~LungMeta,data=CRC_immun)
p.val<-round(Pvalue,3)
HR <- paste("Hazard Ratio= ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low,2), round(up,2), sep = " - "), sep = "")
#if (!dir.exists("metastasis_COHORT")){
#  dir.create("./metastasis_COHORT")
# }
#if(p.val<0.05){
Clinicloutcome<-"Figure1_Clinicaldata"

pdf(file.path(Clinicloutcome,"CRC_Sir_OS——Lung.pdf"),width = 6.49, height = 6.58,onefile = F )
print(ggsurvplot(kmfit,#surv.median.line = "hv", # Add medians survival
                 # Change legends: title & labels
                 main = "Survival curve",
                 #legend.title = paste0(set,"-",sig),
                 legend.labs = c("NoLungMeta","LungMeta"),
                 #xlim = c(0, 2000),
                 xlab="Time(months)",
                 ylab="OS(%)",
                 size = 1,
                 #fun="cumhaz",
                 #fun='event',
                 # Add p-value and tervals
                 pval = TRUE,
                 #test.for.trend = TRUE,###group more than 2groups
                 break.time.by = 6,
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
                 font.tickslab =15,
                 #在左下???标出pvalue、HR???95% CI
                 #???小的p value标为p < 0.001
                 #pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                 # paste("p = ",round(p.val,3), sep = "")), HR, CI, sep = "\n")
) 
)
invisible(dev.off())
###
#####
CRC_immun$metagroup<-ifelse(CRC_immun$Livermeta2==2 & CRC_immun$LungMeta==2,"Liver+Lung",
                            ifelse(CRC_immun$Livermeta2==2 & CRC_immun$LungMeta==1, "Liver",
                                   ifelse(CRC_immun$Livermeta2==1 & CRC_immun$LungMeta==2, "Lung","Others")))

table(CRC_immun$metagroup,CRC_immun$PRSDPD)
barplot(prop.table(table(CRC_immun$疗效评估A1PRSD2PD,CRC_immun$metagroup)),col =c("black","#E31A1C") )

responseTable<-prop.table(table(CRC_immun$疗效评估A1PRSD2PD,CRC_immun$metagroup)) %>% data.frame()
colnames(responseTable)

CRC_immun$response<-ifelse(CRC_immun$疗效评估A1PRSD2PD==1,1,0)
CRC_immun$response<-factor(CRC_immun$response)
ggplot(CRC_immun, aes(x =metagroup, fill=response)) +
  geom_bar(position = "fill",width = 0.4) +
  #scale_fill_manual(name="cat_totalscore",values = c("High"="#EE0000CC","Low"="#008B45CC")) + 
  scale_y_continuous(labels = percent) +
  theme_classic()+
  ylab("response rate")+
  scale_fill_manual(values = c("black","#E31A1C"))+
  geom_text(aes( y=..count../tapply(..count.., ..x.. ,sum)[..x..], 
                 label=paste(..count..,paste("(",percent(..count../tapply(..count.., ..x.. ,sum)[..x..]),")",sep=""),sep=" ")),
            stat="count",color="white", position=position_fill(0.5), vjust=0.5)+
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

#save
ggsave(file.path(Clinicloutcome,"CRC_Sir_response.pdf"),width = 10,height =5)
###cox regression.

rt<-CRC_immun
rt<-data.frame(rt)
colnames(rt)
 vars<-c("ALC","PAB", "LDH","CRP", "CA199","CEA","CA242" )
for(var in vars){
  rt[,var]<-as.numeric(rt[,var])
  rt[,var]<-ifelse(rt[,var]>median(rt[,var]),1,0)
}

Basurv<-Surv(time= rt$OS,event=rt$OS.event)
Univcox<-function(x){
  fml<-as.formula(paste0('Basurv~',x))
  Gcox<-coxph(fml,data=rt)
  GSum<-summary(Gcox)
  HR<-round(GSum$coefficients[,2],3)
  Pvalue<-round(GSum$coefficients[,5],3)
  CI<-paste0(round(GSum$conf.int[,3:4],3),collapse='-')
  Unicox<-data.frame('Characteristics'= rownames(GSum$conf.int),
                     'Characteristics1'= x,
                     'Hazard Ratio'= HR,
                     'CI95'=CI,
                     'Pvalue'=Pvalue)
  return(Unicox)
}
Univcox("Age")
colnames(rt)
rt<-data.frame(rt)
###
##
Varnames<-colnames(rt)[c(13:22,23,31:49)]#
Univar<-lapply(Varnames,Univcox)
Univar<-ldply(Univar,data.frame)

write.csv(Univar,file=file.path(Clinicloutcome,"univarforLM.csv"))

Univar<-Univar[c(1:8,11:16),]

barplot(Univar$Hazard.Ratio)
p <- ggplot(Univar,aes(Characteristics1,Hazard.Ratio))+ coord_flip()
p+geom_bar(stat = 'identity')

sortdf<-Univar[order(Univar$Hazard.Ratio,decreasing = F),]
sortdf$group<-ifelse(sortdf$Pvalue<0.05,1,0)
sortdf$group<-factor(sortdf$group)
sortdf$ID<-sortdf$Characteristics1
sortdf$score<-sortdf$Hazard.Ratio-1
sortdf$score<-as.numeric(sortdf$score)
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)

ggplot(sortdf, aes(ID, score,fill=group)) + 
  geom_bar(stat = 'identity') + 
  coord_flip() + 
  scale_fill_manual(values = c("grey","#E31A1C"), guide = FALSE) + 
  scale_colour_manual(values = c("grey","#E31A1C"), guide = FALSE) +
  #画2条虚线
  geom_hline(yintercept =0, 
            color="black",
             linetype = 2, #画虚线
             size = 0.3) + #线的粗细
  
  #写label
  #geom_text(data = sortdf,
      #     aes(x=ID, y=-0.2, label= paste0(" ", ID), color = group),#bar跟坐标轴间留出间隙
        #    size = 4, #字的大小
         #   hjust = "outward" ) +  #字的对齐方式
  
  xlab("Variables") +ylab("Variable Importance(HR-1)")+
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_rect(size = 1.8)) + #边框粗细
  theme(axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        #axis.text.y = element_blank()
        ) #去除y轴

ggsave(file.path(Clinicloutcome,"CRCvariableImportance.pdf"),width=3.6,height=3.5)

univars<-Univar$Characteristics[Univar$Pvalue<0.05]
###mult cox#
###
Varname2<-colnames(rt)[c(15,50)]#
multfml<-as.formula(paste0('Basurv~',paste0(Varname2,collapse='+')))
#multfml<-as.formula(paste0('Basurv~',paste0(Univar$Characteristics1[Univar$Pvalue<0.05],collapse='+')))
multcox<-coxph(multfml,data=rt)
multisum<-summary(multcox)
multiname<-as.character(univars)
multiHR<-round(multisum$coefficients[,2],3)
multiPvalues<-round(multisum$coefficients[,5],3)
multicilow<-round(multisum$conf.int[,3],3)
multiciup<-round(multisum$conf.int[,4],3)
multi95CI<-paste0(multicilow,'-',multiciup)
mulcox<-data.frame(#'Characteristics'=multiname,
  'Hazard ratio'=multiHR,
  'CI95'=multi95CI,
  'Pvalues'=multiPvalues)

rownames(Univar)<-Univar$Characteristics
mulcox$Characteristics<-rownames(mulcox)
UniMulti<-merge(Univar,mulcox,by="Characteristics",all=TRUE)

write.table(UniMulti,file="UniMulti.txt",sep="\t",quote=F,col.names=T) 

###########非 免疫治疗人群的

CRC_no_immun<- read_excel("Figure1_Clinicaldata/In_house_immunotheapy.xlsx",2)

tmp<-CRC_no_immun

Gcox1<-coxph(Surv(OS,OS.event)~Livermeta,data=tmp)
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
  'P-value'=Pvalue
)
cox.p

##plot the km curves
kmfit<- survfit(Surv(OS,OS.event)~Livermeta,data=tmp)
p.val<-round(Pvalue,3)
HR <- paste("Hazard Ratio= ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low,2), round(up,2), sep = " - "), sep = "")
#if (!dir.exists("metastasis_COHORT")){
#  dir.create("./metastasis_COHORT")
# }


pdf(file.path(Clinicloutcome,"CRC_Sir_OS_no_immune.pdf"),width = 6.49, height = 6.58,onefile = F )
print(ggsurvplot(kmfit,#surv.median.line = "hv", # Add medians survival
                 # Change legends: title & labels
                 main = "Survival curve",
                 #legend.title = paste0(set,"-",sig),
                 legend.labs = c("NoLiverMeta","LiverMeta"),
                 #xlim = c(0, 2000),
                 xlab="Time(months)",
                 ylab="OS(%)",
                 size = 1,
                 #fun="cumhaz",
                 #fun='event',
                 # Add p-value and tervals
                 pval = TRUE,
                 #test.for.trend = TRUE,###group more than 2groups
                 break.time.by = 6,
                 #conf.int = TRUE,
                 #group.by=,
                 # Add risk table
                 risk.table = TRUE,
                 tables.height = 0.185,
                 tables.theme = theme_cleantable(),
                 palette = c("black","darkblue"),
                 #ggtheme = theme_bw(), # Change ggplot2 theme
                 #font.title="OS",
                 font.main =15,
                 font.x =  15,
                 font.y = 15,
                 font.tickslab =15,
                 #在左下???标出pvalue、HR???95% CI
                 #???小的p value标为p < 0.001
                 #pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                 # paste("p = ",round(p.val,3), sep = "")), HR, CI, sep = "\n")
) 
)
invisible(dev.off())
#### for table 


data<-CRC_immun
colnames(data)



Varname<-colnames(data)[c(4,6,7,8,9,10,11,12,13,14,15,16,17,18,20,22,23,24,28,29:48)]
varsToFactor <- colnames(data)[c(4,6,8,10,12,13,16,17,18,20,22,23,24,28)]
data[varsToFactor] <- lapply(data[varsToFactor], factor)

tab<- CreateTableOne(vars = Varname, strata =c("Livermeta"), data = data, test = TRUE,testExact = fisher.test,addOverall = TRUE,
                     testApprox = chisq.test)
tab

nonnormal<-colnames(data)[c(29:48)]

unmatched<-print(tab, 
                                  nonnormal =nonnormal, 
                 showAllLevels = TRUE,smd = TRUE,quote=TRUE)

write.table(unmatched,file="~/Projects/ImmunePro/Figure1/CRC_table_baseline.txt",sep="\t",quote=F,col.names=T) 

