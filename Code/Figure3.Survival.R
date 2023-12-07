##Figure3.Sur
source("/users/zhuxl/HAPS/Code_Script/0.suppressMessages.R")
setwd("/users/zhuxl/HAPS/Code_Script")

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

data$group<-rep("",dim(data)[1])
data[which(data$HAPS>10),]$group<-"High HAPS"
data[which(data$HAPS<=10),]$group<-"Low HAPS"

########################
data$Group1<-paste0(data$group," & ",data$LOH)
data$Group2<-data$Group1
data[which(data$Group2!="High HAPS & LOH=0"),]$Group2<-"Low HAPS | LOH=1"
table(data$Group1);table(data$Group2)

data1<-data[which(data$LOH=="LOH=0"),]
my.surv <- Surv(data1$OS_Months, data1$OS_Event)
fit <- survfit(my.surv ~ data1$group)
data1.survdiff <- survdiff(my.surv ~ data1$group) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))

WES_OS_0<-ggsurvplot(fit, data =data1 ,
           ggtheme=custom_theme(),
           conf.int = F,
           conf.int.style = "step",
           censor = T,size=1.2,censor.size=8,
           title="ICI WES HLA-intact",
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
           paste("p = ",round(p.val,4), sep = ""),
           paste("HR = ",round(HR,2),sep = ""), 
           paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = ""),
           sep = "\n"),
           pval.coord=c(-2, 0.15))

data1<-data[which(data$LOH=="LOH=1"),]
my.surv <- Surv(data1$OS_Months, data1$OS_Event)
fit <- survfit(my.surv ~ data1$group)
data1.survdiff <- survdiff(my.surv ~ data1$group) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))

WES_OS_1<-ggsurvplot(fit, data =data1 ,
           ggtheme=custom_theme(),
           conf.int = F,
           conf.int.style = "step",
           censor = T,size=1.2,censor.size=8,
           title="ICI WES HLA-LOH",
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
           paste("p = ",round(p.val,4), sep = ""),
           paste("HR = ",round(HR,2),sep = ""), 
           paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = ""),
           sep = "\n"),
           pval.coord=c(-2, 0.15))

data1<-data
my.surv <- Surv(data1$OS_Months, data1$OS_Event)
fit <- survfit(my.surv ~ data1$Group1)
data1.survdiff <- survdiff(my.surv ~ data1$Group1) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))


WES_OS_Group1<-ggsurvplot(fit, data =data1 ,
           ggtheme=custom_theme(),
           conf.int = F,
           conf.int.style = "step",
           censor = T,size=1.2,censor.size=8,
           title="ICI WES HAPS & HLA",
           ylab="Overall Survival (%)",
           xlab="Time (Months)",
           palette = "npg", 
           risk.table=TRUE, 
           risk.table.title="",
           risk.table.y.text=T,
           risk.table.height=0.3,
           font.legend = 13,       
           legend.title="",
           legend="none",
           legend.labs=c("High HAPS & HLA-intact",
           "High HAPS & HLA-LOH",
           "Low HAPS & HLA-intact",
           "Low HAPS & HLA-LOH"),
           font.x = 0,font.y = 13,font.tickslab=13,
           pval = paste(
           paste("p = ",round(p.val,4), sep = ""),
           paste("HR = ",round(HR,2),sep = ""), 
           paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = ""),
           sep = "\n"),
           pval.coord=c(-2, 0.15))



data1<-data
my.surv <- Surv(data1$OS_Months, data1$OS_Event)
fit <- survfit(my.surv ~ data1$Group2)
data1.survdiff <- survdiff(my.surv ~ data1$Group2) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))

WES_OS_Group2<-ggsurvplot(fit, data =data1 ,
           ggtheme=custom_theme(),
           conf.int = F,
           conf.int.style = "step",
           censor = T,size=1.2,censor.size=8,
           title="ICI WES HAPS & HLA",
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
           legend.labs=c("High HAPS & HLA-intact",
           "Low HAPS | HLA-LOH"),
           font.x = 0,font.y = 13,font.tickslab=13,
           pval = paste(
           paste("p = ",round(p.val,4), sep = ""),
           paste("HR = ",round(HR,2),sep = ""), 
           paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = ""),
           sep = "\n"),
           pval.coord=c(-2, 0.15))


##################################
data<-read.xlsx("./WES_Cohort",sheet=1)
data<-data[which(data$PFS_Event==0 |
data$PFS_Event==1),]

data$Group1<-paste0(data$group," & ",data$LOH)
data$Group2<-data$Group1
data[which(data$Group2!="High HAPS & LOH=0"),]$Group2<-"Low HAPS | LOH=1"
table(data$Group1);table(data$Group2)

data1<-data[which(data$LOH=="LOH=0"),]
my.surv <- Surv(data1$PFS_Months, data1$PFS_Event)
fit <- survfit(my.surv ~ data1$group)
data1.survdiff <- survdiff(my.surv ~ data1$group) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))

WES_PFS_0<-ggsurvplot(fit, data =data1 ,
           ggtheme=custom_theme(),
           conf.int = F,
           conf.int.style = "step",
           censor = T,size=1.2,censor.size=8,
           title="ICI WES HLA-intact",
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
           paste("p = ",round(p.val,4), sep = ""),
           paste("HR = ",round(HR,2),sep = ""), 
           paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = ""),
           sep = "\n"),
           pval.coord=c(-2, 0.15))


data1<-data[which(data$LOH=="LOH=1"),]
my.surv <- Surv(data1$PFS_Months, data1$PFS_Event)
fit <- survfit(my.surv ~ data1$group)
data1.survdiff <- survdiff(my.surv ~ data1$group) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))

WES_PFS_1<-ggsurvplot(fit, data =data1 ,
           ggtheme=custom_theme(),
           conf.int = F,
           conf.int.style = "step",
           censor = T,size=1.2,censor.size=8,
           title="ICI WES HLA-LOH",
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
           paste("p = ",round(p.val,4), sep = ""),
           paste("HR = ",round(HR,2),sep = ""), 
           paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = ""),
           sep = "\n"),
           pval.coord=c(-2, 0.15))


data1<-data
my.surv <- Surv(data1$PFS_Months, data1$PFS_Event)
fit <- survfit(my.surv ~ data1$Group1)
data1.survdiff <- survdiff(my.surv ~ data1$Group1) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))

WES_PFS_Group1<-ggsurvplot(fit, data =data1 ,
           ggtheme=custom_theme(),
           conf.int = F,
           conf.int.style = "step",
           censor = T,size=1.2,censor.size=8,
           title="ICI WES HAPS & HLA",
           ylab="Progression Free Survival (%)",
           xlab="Time (Months)",
           palette = "npg", 
           risk.table=TRUE, 
           risk.table.title="",
           risk.table.y.text=T,
           risk.table.height=0.3,
           font.legend = 13,       
           legend.title="",
           legend="none",
           legend.labs=c("High HAPS & HLA-intact",
           "High HAPS & HLA-LOH",
           "Low HAPS & HLA-intact",
           "Low HAPS & HLA-LOH"),
           font.x = 0,font.y = 13,font.tickslab=13,
           pval = paste(
           paste("p = ",round(p.val,4), sep = ""),
           paste("HR = ",round(HR,2),sep = ""), 
           paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = ""),
           sep = "\n"),
           pval.coord=c(-2, 0.15))



data1<-data
my.surv <- Surv(data1$PFS_Months, data1$PFS_Event)
fit <- survfit(my.surv ~ data1$Group2)
data1.survdiff <- survdiff(my.surv ~ data1$Group2) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))

WES_PFS_Group2<-ggsurvplot(fit, data =data1 ,
           ggtheme=custom_theme(),
           conf.int = F,
           conf.int.style = "step",
           censor = T,size=1.2,censor.size=8,
           title="ICI WES HAPS & HLA",
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
           legend.labs=c("High HAPS & HLA-intact",
           "Low HAPS | HLA-LOH"),
           font.x = 0,font.y = 13,font.tickslab=13,
           pval = paste(
           paste("p = ",round(p.val,4), sep = ""),
           paste("HR = ",round(HR,2),sep = ""), 
           paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = ""),
           sep = "\n"),
           pval.coord=c(-2, 0.15))

####################################

##########################
data<-read.xlsx("./Wang_Panel.xlsx",sheet=1)
data1<-data[grep("Wang-Panel-T",data$Cohort),]

data1<-data[which(data$LOH=="LOH=0"),]
my.surv <- Surv(data1$OS_Months, data1$OS_Event)
fit <- survfit(my.surv ~ data1$Group)
data1.survdiff <- survdiff(my.surv ~ data1$Group) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))

Panel_T_0<-ggsurvplot(fit, data =data1 ,
           gtheme=custom_theme(),
           conf.int = F,
           conf.int.style = "step",
           censor = T,size=1.2,censor.size=8,
           title="Wang-Panel-T HLA-intact",
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
           paste("p = ",round(p.val,4), sep = ""),
           paste("HR = ",round(HR,2),sep = ""), 
           paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = ""),
           sep = "\n"),
           pval.coord=c(-2, 0.15))

data1<-data[which(data$LOH=="LOH=1"),]
my.surv <- Surv(data1$OS_Months, data1$OS_Event)
fit <- survfit(my.surv ~ data1$Group)
data1.survdiff <- survdiff(my.surv ~ data1$Group) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))

Panel_T_1<-ggsurvplot(fit, data =data1 ,
           ggtheme=custom_theme(),
           conf.int = F,
           conf.int.style = "step",
           censor = T,size=1.2,censor.size=8,
           title="Wang-Panel-T HLA-LOH",
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
           paste("p = ",round(p.val,4), sep = ""),
           paste("HR = ",round(HR,2),sep = ""), 
           paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = ""),
           sep = "\n"),
           pval.coord=c(-2, 0.15))

data1<-data
my.surv <- Surv(data1$OS_Months, data1$OS_Event)
fit <- survfit(my.surv ~ data1$HAPS_LOH_4)
data1.survdiff <- survdiff(my.surv ~ data1$HAPS_LOH_4) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))

Panel_T_Group1<-ggsurvplot(fit, data =data1 ,
           ggtheme=custom_theme(),
           conf.int = F,
           conf.int.style = "step",
           censor = T,size=1.2,censor.size=8,
           title="Wang-Panel-T HAPS & HLA",
           ylab="Overall Survival (%)",
           xlab="Time (Months)",
           palette = "npg", 
           risk.table=TRUE, 
           risk.table.title="",
           risk.table.y.text=T,
           risk.table.height=0.3,
           font.legend = 13,       
           legend.title="",
           legend="none",
           legend.labs=c("High HAPS & HLA-intact",
           "High HAPS & HLA-LOH",
           "Low HAPS & HLA-intact",
           "Low HAPS & HLA-LOH"),
           font.x = 0,font.y = 13,font.tickslab=13,
           pval = paste(
           paste("p = ",round(p.val,4), sep = ""),
           paste("HR = ",round(HR,2),sep = ""), 
           paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = ""),
           sep = "\n"),
           pval.coord=c(-2, 0.15))


data1<-data
my.surv <- Surv(data1$OS_Months, data1$OS_Event)
fit <- survfit(my.surv ~ data1$HAPS_LOH_2)
data1.survdiff <- survdiff(my.surv ~ data1$HAPS_LOH_2) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))

Panel_T_Group2<-ggsurvplot(fit, data =data1 ,
           ggtheme=custom_theme(),
           conf.int = F,
           conf.int.style = "step",
           censor = T,size=1.2,censor.size=8,
           title="Wang-Panel-T HAPS & HLA",
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
           legend.labs=c("High HAPS & HLA-intact",
           "Low HAPS | HLA-LOH"),
           font.x = 0,font.y = 13,font.tickslab=13,
           pval = paste(
           paste("p = ",round(p.val,4), sep = ""),
           paste("HR = ",round(HR,2),sep = ""), 
           paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = ""),
           sep = "\n"),
           pval.coord=c(-2, 0.15))

##########################
data<-read.xlsx("./Wang_Panel.xlsx",sheet=1)
data1<-data[grep("Wang-Panel-B",data$Cohort),]

data1<-data[which(data$LOH=="LOH=0"),]
my.surv <- Surv(data1$OS_Months, data1$OS_Event)
fit <- survfit(my.surv ~ data1$Group)
data1.survdiff <- survdiff(my.surv ~ data1$Group) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))

Panel_B_0<-ggsurvplot(fit, data =data1 ,
           #surv.median.line = "hv",
           #ggtheme = theme(plot.title=element_text(hjust=0.5)),
           ggtheme=custom_theme(),
           conf.int = F,
           conf.int.style = "step",
           censor = T,size=1.2,censor.size=8,
           title="Wang-Panel-B HLA-intact",
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
           paste("p = ",round(p.val,4), sep = ""),
           paste("HR = ",round(HR,2),sep = ""), 
           paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = ""),
           sep = "\n"),
           pval.coord=c(-2, 0.15))

data1<-data[which(data$LOH=="LOH=1"),]
my.surv <- Surv(data1$OS_Months, data1$OS_Event)
fit <- survfit(my.surv ~ data1$Blood)
data1.survdiff <- survdiff(my.surv ~ data1$Blood) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))

Panel_B_1<-ggsurvplot(fit, data =data1 ,
           #surv.median.line = "hv",
           #ggtheme = theme(plot.title=element_text(hjust=0.5)),
           ggtheme=custom_theme(),
           conf.int = F,
           conf.int.style = "step",
           censor = T,size=1.2,censor.size=8,
           title="Wang-Panel-B HLA-LOH",
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
           paste("p = ",round(p.val,4), sep = ""),
           paste("HR = ",round(HR,2),sep = ""), 
           paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = ""),
           sep = "\n"),
           pval.coord=c(-2, 0.15))

data1<-data
my.surv <- Surv(data1$OS_Months, data1$OS_Event)
fit <- survfit(my.surv ~ data1$B_L)
data1.survdiff <- survdiff(my.surv ~ data1$B_L) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))

Panel_B_Group1<-ggsurvplot(fit, data =data1 ,
           #surv.median.line = "hv",
           #ggtheme = theme(plot.title=element_text(hjust=0.5)),
           ggtheme=custom_theme(),
           conf.int = F,
           conf.int.style = "step",
           censor = T,size=1.2,censor.size=8,
           title="Wang-Panel-B HAPS & HLA",
           ylab="Overall Survival (%)",
           xlab="Time (Months)",
           palette = "npg", 
           risk.table=TRUE, 
           risk.table.title="",
           risk.table.y.text=T,
           risk.table.height=0.3,
           font.legend = 13,       
           legend.title="",
           legend="none",
           legend.labs=c("High HAPS & HLA-intact",
           "High HAPS & HLA-LOH",
           "Low HAPS & HLA-intact",
           "Low HAPS & HLA-LOH"),
           font.x = 0,font.y = 13,font.tickslab=13,
           pval = paste(
           paste("p = ",round(p.val,4), sep = ""),
           paste("HR = ",round(HR,2),sep = ""), 
           paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = ""),
           sep = "\n"),
           pval.coord=c(-2, 0.15))


data1<-data
my.surv <- Surv(data1$OS_Months, data1$OS_Event)
fit <- survfit(my.surv ~ data1$B_L_2)
data1.survdiff <- survdiff(my.surv ~ data1$B_L_2) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))

Panel_B_Group2<-ggsurvplot(fit, data =data1 ,
           #surv.median.line = "hv",
           #ggtheme = theme(plot.title=element_text(hjust=0.5)),
           ggtheme=custom_theme(),
           conf.int = F,
           conf.int.style = "step",
           censor = T,size=1.2,censor.size=8,
           title="Wang-Panel-B HAPS & HLA",
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
           legend.labs=c("High HAPS & HLA-intact",
           "Low HAPS | HLA-LOH"),
           font.x = 0,font.y = 13,font.tickslab=13,
           pval = paste(
           paste("p = ",round(p.val,4), sep = ""),
           paste("HR = ",round(HR,2),sep = ""), 
           paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = ""),
           sep = "\n"),
           pval.coord=c(-2, 0.15))

##########################

p<-arrange_ggsurvplots(list(
WES_OS_0,WES_PFS_0,Panel_T_0,Panel_B_0,
WES_OS_1,WES_PFS_1,Panel_T_1,Panel_B_1,
WES_OS_Group1,WES_PFS_Group1,Panel_T_Group1,Panel_B_Group1,
WES_OS_Group2,WES_PFS_Group2,Panel_T_Group2,Panel_B_Group2
),
nrow=4,ncol=4)
ggsave("Figure3.WES_Panel_4_Sur.pdf",p,height = 20,width = 25)



