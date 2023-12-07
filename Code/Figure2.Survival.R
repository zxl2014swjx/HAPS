##Figure2.Sur
source("/users/zhuxl/HAPS/Code_Script/0.suppressMessages.R")
setwd("/users/zhuxl/HAPS/Code_Script/")

data<-read.xlsx("./WES_Cohort.xlsx",sheet=1)
dim(data)
#[1] 1125   29
#head(data)
colnames(data)
# [1] "Cohort"       "Batch"        "Sample"      
# [4] "Cancer"       "LOH"          "TMB"         
# [7] "Smoking"      "PD_L1"        "RECIST"      
#[10] "OS_Months"    "OS_Event"     "PFS_Months"  
#[13] "PFS_Event"    "Homozygous"   "HLAI.Alleles"
#[16] "HLA.A"        "HLA.B"        "HLA.C"       
#[19] "HLA"          "TNB.A"        "TNB.B"       
#[22] "TNB.C"        "TNB"          "Age"         
#[25] "Gender"       "HAPS.A"       "HAPS.B"      
#[28] "HAPS.C"       "HAPS"        

custom_theme<-function(){
  theme_survminer()%+replace%
  theme(plot.title=element_text(hjust=0.5))
}

data$group<-rep("",dim(data)[1])
data[which(data$HAPS>10),]$group<-"High HAPS"
data[which(data$HAPS<=10),]$group<-"Low HAPS"

li<-list()
data<-data[-grep("Training",data$Batch),]
var<-as.data.frame(table(data$Batch))
#         Var1 Freq
#1        TCGA  333
#2 Validation1  249
#3 Validation2  479

for(i in 1:dim(var)[1]){
data1<-data[which(data$Batch==var[i,1]),]
my.surv <- Surv(data1$OS_Months, data1$OS_Event)
fit <- survfit(my.surv ~ data1$group)
data1.survdiff <- survdiff(my.surv ~ data1$group) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))

li[[i]]<-ggsurvplot(fit, data =data1 ,
           ggtheme=custom_theme(),
           conf.int = F,
           conf.int.style = "step",
           censor = T,
           size=1.2,censor.size=8,
           title=as.vector(var[i,1]),
           ylab="Overall Survival (%)",
           xlab="Time (Months)",
           palette = c("#CB3425","#3F5688"), 
           risk.table=TRUE, 
           risk.table.title="",
           risk.table.y.text=T,
           risk.table.height=0.3,
           font.legend = 13,       
           legend.title="",
           legend="none",
           legend.labs=c("High HAPS","Low HAPS"),
           font.x = 0,font.y = 13,font.tickslab=13,
           pval = paste(
           paste("P = ",round(p.val,4), sep = ""),
           paste("HR = ",round(HR,2),sep = ""), 
           paste("95%CI: ", paste(round(low95,2), 
           round(up95,2), sep = " - "), sep = ""),
           sep = "\n"),
           pval.coord=c(-2, 0.15))


}

data<-read.xlsx("../WES_Cohort.xlsx",sheet=1)
data$group<-rep("",dim(data)[1])
data[which(data$HAPS>10),]$group<-"High HAPS"
data[which(data$HAPS<=10),]$group<-"Low HAPS"

data2<-data[which(data$PFS_Event==0 |
data$PFS_Event==1),]
dim(data2)
#[1] 488  30
table(data2$PFS_Event)
#  0   1 
#217 271 

my.surv <- Surv(data2$PFS_Months, data2$PFS_Event)
fit <- survfit(my.surv ~ data2$group)
data1.survdiff <- survdiff(my.surv ~ data2$group) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))

PFS<-ggsurvplot(fit, data =data2 ,
           ggtheme=custom_theme(),
           conf.int = F,
           conf.int.style = "step",
           censor = T,
           size=1.2,
           censor.size=8,
           title="PFS",
           ylab="Progression Free Survival (%)",
           xlab="Time (Months)",
           palette = c("#CB3425","#3F5688"), 
           risk.table=TRUE, 
           risk.table.title="",
           risk.table.y.text=T,
           risk.table.height=0.3,
           font.legend = 13,       
           legend.title="",
           legend="none",
           legend.labs=c("High HAPS","Low HAPS"),
           font.x = 0,font.y = 13,font.tickslab=13,
           pval = paste(
           paste("P = ",round(p.val,4), sep = ""),
           paste("HR = ",round(HR,2),sep = ""), 
           paste("95%CI: ", paste(round(low95,2), 
           round(up95,2), sep = " - "), sep = ""),
           sep = "\n"),
           pval.coord=c(-2, 0.15))


############
data<-read.xlsx("./Wang_Panel.xlsx",sheet=1)
dim(data)
#[1] 93 24

data1<-data[grep("Wang-Panel-T",data$Cohort),]
my.surv <- Surv(data1$OS_Months, data1$OS_Event)
fit <- survfit(my.surv ~ data1$Group)
data1.survdiff <- survdiff(my.surv ~ data1$Group) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))

Tissue<-ggsurvplot(fit, data =data1 ,
           ggtheme=custom_theme(),
           conf.int = F,
           conf.int.style = "step",
           censor = T,
           size=1.2,censor.size=8,
           title="Wang-Panel-T cohort",
           ylab="Overall Survival (%)",
           xlab="Time (Months)",
           palette = c("#CB3425","#3F5688"), 
           risk.table=TRUE, 
           risk.table.title="",
           risk.table.y.text=T,
           risk.table.height=0.3,
           font.legend = 13,       
           legend.title="",
           legend="none",
           legend.labs=c("High HAPS","Low HAPS"),
           font.x = 0,font.y = 13,font.tickslab=13,
           pval = paste(
           paste("P = ",round(p.val,3), sep = ""),
           paste("HR = ",round(HR,2),sep = ""), 
           paste("95%CI: ", paste(round(low95,2), 
           round(up95,2), sep = " - "), sep = ""),
           sep = "\n"),
           pval.coord=c(-2, 0.15))


data1<-data[grep("Wang-Panel-B",data$Cohort),]
my.surv <- Surv(data1$OS_Months, data1$OS_Event)
fit <- survfit(my.surv ~ data1$Group)
data1.survdiff <- survdiff(my.surv ~ data1$Group) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))

Blood<-ggsurvplot(fit, data =data1 ,
           ggtheme=custom_theme(),
           conf.int = F,
           conf.int.style = "step",
           censor = T,
           size=1.2,censor.size=8,
           title="Wang-Panel-B cohort",
           ylab="Overall Survival (%)",
           xlab="Time (Months)",
           palette = c("#CB3425","#3F5688"), 
           risk.table=TRUE, 
           risk.table.title="",
           risk.table.y.text=T,
           risk.table.height=0.3,
           font.legend = 13,       
           legend.title="",
           legend="none",
           legend.labs=c("High HAPS","Low HAPS"),
           font.x = 0,font.y = 13,font.tickslab=13,
           pval = paste(
           paste("P = ",round(p.val,3), sep = ""),
           paste("HR = ",round(HR,2),sep = ""), 
           paste("95%CI: ", paste(round(low95,2), 
           round(up95,2), sep = " - "), sep = ""),
           sep = "\n"),
           pval.coord=c(-2, 0.15))


p<-arrange_ggsurvplots(list(li[[2]],PFS,Tissue,
li[[3]],li[[1]],Blood),
nrow=3,ncol=2)
dev.off()
dev.off()
ggsave("Figure2.Survival.pdf",p,height = 12,width = 10)

###############
