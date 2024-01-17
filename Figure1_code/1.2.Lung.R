###
library(plyr)
library(survival)
library(readxl)
Lung_immun<- read_excel("~/Projects/ImmunePro/Figure1/In_house_immunotheapy.xlsx",3)
#
Lung_immun$PFS<-as.numeric(Lung_immun$PFS)
tmp<-Lung_immun[Lung_immun$surgery==0,]
citPie <- prop.table(table(tmp$Livermeta))
print(citPie)
pie(citPie,
    col=c("grey","#E31A1C"))

citPie <- prop.table(table(tmp$`Brain metastasis`))
print(citPie)
pie(citPie,
    col=c("grey","#E31A1C"))

citPie <- prop.table(table(tmp$`Bone transfer`))
print(citPie)
pie(citPie,
    col=c("grey","#E31A1C"))

citPie <- prop.table(table(tmp$`adrenal gland`))
print(citPie)
pie(citPie,
    col=c("grey","#E31A1C"))

tmp<-data.frame(tmp)
citpies<-c()
vars<-c("Livermeta","Brain.metastasis","adrenal.gland","Bone.transfer")

for(var in vars){
  citPie <- prop.table(table(tmp[,var]))
  citpies<-rbind(citpies,citPie)
}
##
rownames(citpies)<-vars
sotmpdf<-citpies
sotmpdf<-data.frame(sotmpdf)
sotmpdf$ID<-rownames(sotmpdf)
sotmpdf$score<-sotmpdf$Yes


ggplot(sotmpdf, aes(ID, score,fill=ID)) + 
  geom_bar(stat = 'identity') + 
  coord_flip() + 
  scale_fill_manual(values = c("grey","grey","grey","#E31A1C","grey"), guide = FALSE) + 
  scale_colour_manual(values =c("grey","grey","grey","#E31A1C","grey"), guide = FALSE) +
  #画2条虚线
  geom_hline(yintercept =0.2, 
             color="black",
             linetype = 2, #画虚线
             size = 0.3) + #线的粗细
  
  #写label
  #geom_text(data = sotmpdf,
  #     aes(x=ID, y=-0.2, label= paste0(" ", ID), color = group),#bar跟坐标轴间留出间隙
  #    size = 4, #字的大小
  #   hjust = "outward" ) +  #字的对齐方式
  
  xlab("Metastasis Sits") +ylab("% Propotmpion")+
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_rect(size = 1.8)) + #边框粗细
  theme(axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        #axis.text.y = element_blank()
  ) #去除y轴
###
ggsave(file.path(Clinicloutcome,"Lung_metavariableImportance.pdf"),width=3.6,height=3.5)

tmp$number.of.transfers<-ifelse(tmp$number.of.transfers>1,1,0)

Basurv<-Surv(time= tmp$OS,event=tmp$OS.event)
Univcox<-function(x){
  fml<-as.formula(paste0('Basurv~',x))
  Gcox<-coxph(fml,data=tmp)
  GSum<-summary(Gcox)
  HR<-round(GSum$coefficients[,2],3)
  Pvalue<-round(GSum$coefficients[,5],3)
  CI<-paste0(round(GSum$conf.int[,3:4],3),collapse='-')
  Unicox<-data.frame('Characteristics'= rownames(GSum$conf.int),
                     'Characteristics1'= x,
                     'HR'= HR,
                     'CI95'=CI,
                     "low"=round(GSum$conf.int[,3],3),
                     "up"=round(GSum$conf.int[,4],3),
                     'Pvalue'=Pvalue)
  return(Unicox)
}
Univcox("age")
colnames(tmp)
tmp<-data.frame(tmp)
###
Varnames<-c("number.of.transfers","Livermeta","Brain.metastasis","Bone.transfer","adrenal.gland")
Varnames<-colnames(tmp)[c(27:30,31)]#
Univar<-lapply(Varnames,Univcox)
Univar<-ldply(Univar,data.frame)
Univar
#########################forest plot for lung
library(forestplot)
library(readxl)
rt <-Univar
colnames(rt)
rownames(rt)<-rt$Characteristics1
rt$HR<-round(as.numeric(rt$HR),2)
rt$up<-round(as.numeric(rt$up),2)
rt$low<-round(as.numeric(rt$low),2)
rt$pvalue<-round(as.numeric(rt$Pvalue),3)
tabletext <- cbind(c("\nDataset",NA,rownames(rt),NA),
                   c("Hazard Ratio\n(95% CI)",NA, 
                     paste0(format(rt$HR,nsmall=2),
                            " (",format(rt$low,nsmall = 2),"-",format(rt$up,nsmall = 2),")",sep=""),NA),
                   c("p-value",NA,rt$pvalue,NA))
pdf(file.path(Clinicloutcome,paste0("Lung_OS_all","forestplot.pdf")),width = 8,height =8,onefile=FALSE)
#pdf(paste0(cellnames[3],"forestplot_OS.pdf"),width = 8,height = 5,onefile=FALSE)
forestplot(labeltext=tabletext, #
           mean=c(NA,NA,rt$HR,NA),#HR
           lower=c(NA,NA,rt$low,NA), #95%
           upper=c(NA,NA,rt$up,NA),#95%
           #title="Hazard Ratio",
           graph.pos=2,#????????
           graphwidth = unit(.3,"npc"),#
           #fn.ci_norm="fpDrawDiamondCI",#
           col=fpColors(box="steelblue", lines="black", zero = "black"),#
           #col=fpColors(box="black", lines="black", zero = "black"),#
           #boxsize=c(NA,NA,NA,rt$numberofpatients,NA)/200,#
           boxsize=0.5,
           #lwd.ci=2,ci.vertices.height = 0.1,ci.vertices=TRUE,#
           zero=1,#zero
           #xlog = TRUE,
           lwd.zero=1,#zero
           #grid = structure(c(rt[1,]$Hazard.Ratio), gp = gpar(col = "black", lty=2,lwd=2)),#()
           #xticks = c(0.5,0.75, 1,1.25,1.5),#
           #clip = c(0.1,2.5), 
           #lwd.xaxis=2,#X
           #xlab="     <-Favour Combination  Therapy       Favour Target  Therapy->",#X
           hrzl_lines=list("3" = gpar(lwd=2, col="black"),#
                           #"4" = gpar(lwd=60,lineend="butt", columns=c(1:4), col="#99999922"),#
                           "9" = gpar(lwd=2, col="black")),#"nrow(rt)+5
           txt_gp=fpTxtGp(label=gpar(cex=1.25),#
                          ticks=gpar(cex=1.25),
                          xlab=gpar(cex = 1.25),
                          title=gpar(cex = 1.25)),
           #is.summary = c(T,rep(F,27)),#
           #lineheight = unit(.75,"cm"),#
           align=c("l","c","c"),#
           #cex=10,
           colgap = unit(0.1,"cm"),#
           #mar=unit(rep(0.25, times = 4), "cm"),#
           new_page = T#
)
dev.off()
####
colnames(tmp)
##plot the km curves
kmfit<- survfit(Surv(OS,OS.event)~Livermeta+Bone.transfer,data=tmp)
p.val<-round(Pvalue,3)
HR <- paste("Hazard Ratio= ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low,2), round(up,2), sep = " - "), sep = "")
Clinicloutcome<-"~/Projects/ImmunePro/Figure1"
pdf(file.path(Clinicloutcome,"Lung_Sir_OS_bone_liver.pdf"),width = 6.49, height = 6.58,onefile = F )
print(ggsurvplot(kmfit,#surv.median.line = "hv", # Add medians survival
                 # Change legends: title & labels
                 main = "Survival curve",
                 #legend.title = paste0(set,"-",sig),
                 #legend.labs = c("NoLiverMeta","LiverMeta"),
                 legend.labs = c("NoLiver,NoBone","NoLiver,Bone","Liver,NoBone","Liver,Bone"),
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
                 palette = c("black","grey","darkblue","#E31A1C"),
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




data<-tmp
colnames(data)



Varname<-colnames(data)[c(4:9,11,12,13,14,15,16,19,25,41:62,92,93)]
varsToFactor <- colnames(data)[c(5,7,8,9,11,12,13,14)]
data[varsToFactor] <- lapply(data[varsToFactor], factor)

varsTonumber <- colnames(data)[c(15,19,41:62,92,93)]
data[varsTonumber ] <- lapply(data[varsTonumber], as.numeric)


?CreateTableOne
tab<- CreateTableOne(vars = Varname, strata =c("Livermeta"), data = data, test = TRUE,testExact = fisher.test, addOverall = TRUE,smd = TRUE,
                     testApprox = chisq.test)
tab

nonnormal<-colnames(data)[c(41:62,92,93)]


unmatched<-print(tab, 
                 nonnormal =nonnormal, 
                 showAllLevels = TRUE,smd = TRUE,quote=TRUE)

write.table(unmatched,file="~/Projects/ImmunePro/Figure1/Lung_table_baseline.txt",sep="\t",quote=F,col.names=T) 
