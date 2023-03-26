##Figure1.A & C
source("/users/zhuxl/HAPS/Code_Script/0.suppressMessages.R")
setwd("/users/zhuxl/HAPS/Code_Script/Figure 1")

data<-read.xlsx("../WES_Cohort.xlsx",sheet=1)
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

#HLA.ABC(16,17,18)
HLA<-melt(data,id.vars=c(1:15,19:29))
table(HLA$variable)
#HLA.A HLA.B HLA.C 
# 1125  1125  1125 
my_comparisons <- list( c("HLA.B", "HLA.C"),
c("HLA.A", "HLA.B"),
c("HLA.A", "HLA.C"))

HLA_gg<-ggplot(HLA,aes(x = variable, y = value,
fill = variable)) +
geom_violin(position = position_dodge(0.9),
alpha = 0.8,width = 1.2,trim = F,color = "grey50") +
geom_boxplot(width = 0.3,
show.legend = F,position = position_dodge(0.9),
color = 'grey20',alpha = 0.8,
outlier.color = 'grey50') +
#geom_point(color="black",size=1,position="identity")+
#geom_jitter(position="jitter",colour="grey20")+
xlab("")+ylab("HLA-I alleles divergence")+
theme_classic()+
scale_fill_manual(values = c("#E64B35","#4DBBD5","grey50"))+
theme(legend.position="none",
text=element_text(size=15,hjust = 0.5))+
stat_compare_means(comparisons = my_comparisons)

HLA_den<-ggplot(HLA,aes(x=value,fill=variable))+
geom_density(alpha=0.8)+
theme_classic()+
scale_fill_manual(values = c("#E64B35","#4DBBD5","grey50"))+
xlab("HLA-I Alleles Divergence")+
ylab("The Ratio of HLA Divergence")+
theme(legend.title=element_blank(),
text=element_text(size=15,hjust = 0.5),
legend.position="top")


#TNB.ABC(20,21,22)
TNB<-melt(data,id.vars=c(1:19,23:29))
table(TNB$variable)
#TNB.A TNB.B TNB.C 
# 1125  1125  1125 
TNB$value1<-log10(TNB$value+1)
my_comparisons <- list( c("TNB.B", "TNB.C"),
c("TNB.A", "TNB.B"),
c("TNB.A", "TNB.C"))

TNB_gg<-ggplot(TNB,aes(x = variable, y = value1,
fill = variable)) +
geom_violin(position = position_dodge(0.9),
alpha = 0.8,width = 1.2,trim = F,color = "grey50") +
geom_boxplot(width = 0.3,
show.legend = F,position = position_dodge(0.9),
color = 'grey20',alpha = 0.8,
outlier.color = 'grey50') +
#geom_point(color="black",size=1,position="identity")+
#geom_jitter(position="jitter",colour="grey20")+
xlab("")+ylab("log10(Tumor Neoantigen Burden + 1)")+
theme_classic()+
scale_fill_manual(values = c("#E64B35","#4DBBD5","grey50"))+
theme(legend.position="none",
text=element_text(size=15,hjust = 0.5))+
stat_compare_means(comparisons = my_comparisons)

TNB_den<-ggplot(TNB,aes(x=value1,fill=variable))+
geom_density(alpha=0.8)+
theme_classic()+
scale_fill_manual(values = c("#E64B35","#4DBBD5","grey50"))+
xlab("log10(Tumor Neoantigen Burden + 1)")+
ylab("The Ratio of log10(TNB+1)")+
theme(legend.title=element_blank(),
text=element_text(size=15,hjust = 0.5),
legend.position="top")

#HAPS.ABC(26,27,28)
HAPS<-melt(data,id.vars=c(1:25,29))
table(HAPS$variable)
#HAPS.A HAPS.B HAPS.C 
#  1125   1125   1125 
my_comparisons <- list( c("HAPS.B", "HAPS.C"), 
c("HAPS.A", "HAPS.B"),
c("HAPS.A", "HAPS.C"))

HAPS_gg<-ggplot(HAPS,aes(x = variable, y = value,
fill = variable)) +
geom_violin(position = position_dodge(0.9),
alpha = 0.8,width = 1.2,trim = F,color = "grey50") +
geom_boxplot(width = 0.3,
show.legend = F,position = position_dodge(0.9),
color = 'grey20',alpha = 0.8,
outlier.color = 'grey50') +
#geom_point(color="black",size=1,position="identity")+
#geom_jitter(position="jitter",colour="grey20")+
xlab("")+ylab("HLA-I Antigen Presentation Score")+
theme_classic()+
scale_fill_manual(values = c("#E64B35","#4DBBD5","grey50"))+
theme(legend.position="none",
text=element_text(size=15,hjust = 0.5))+
stat_compare_means(comparisons = my_comparisons)

HAPS_den<-ggplot(HAPS,aes(x=value,fill=variable))+
geom_density(alpha=0.8)+
theme_classic()+
scale_fill_manual(values = c("#E64B35","#4DBBD5","grey50"))+
xlab("HLA-I Antigen Presentation Score")+
ylab("The Ratio of HAPS")+
theme(legend.title=element_blank(),
text=element_text(size=15,hjust = 0.5),
legend.position="top")

best1<-plot_grid(HLA_den,TNB_den,HAPS_den,
HLA_gg,TNB_gg,HAPS_gg,
nrow=2,ncol=3,
rel_heights = c(1,1))

ggsave("Figure1.A&C.pdf",best1,height = 9,width =15)


#######################################
#Figure1.D&E up: cutoff
AL<-data.frame()
data<-data[grep("Training",data$Batch),]
var<-as.data.frame(table(data$Batch))
for(i in 1:dim(var)[1]){
Batch<-as.vector(var[i,1])
data1<-data[which(data$Batch==Batch),]
H<-data.frame()
cut<-seq(min(data1$HAPS)+0.1,max(data1$HAPS)-0.1,0.1)
for(j in 1:length(cut)){
cutoff<-cut[j]
data1$group<-rep("",dim(data1)[1])
data1[which(data1$HAPS>cutoff),]$group<-"High HAPS"
data1[which(data1$HAPS<=cutoff),]$group<-"Low HAPS"
my.surv <- Surv(data1$OS_Months, data1$OS_Event)
fit <- survfit(my.surv ~ data1$group)
data1.survdiff <- survdiff(my.surv ~ data1$group) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))
al<-cbind(cutoff,p.val,HR,low95,up95)
H<-rbind(H,al)
}
al<-cbind(Batch,H)
AL<-rbind(AL,al)
}



plt<-melt(AL,id.vars=c(1,2))
li<-list()
var<-as.data.frame(table(AL$Batch))
for(i in 1:dim(var)[1]){
Batch<-as.vector(var[i,1])
AL1<-AL[which(AL$Batch==Batch),]
li[[i]]<-ggplot(AL1,aes(x=cutoff))+
geom_line(aes(y=HR),col="#CB3425",size=1.2)+
geom_line(aes(y=p.val),col="#3F5688",size=1.2)+
ggtitle(Batch)+
theme_classic()+
theme(plot.title = element_text(hjust = 0.5))+
ylab("Red:HR; Blue:Pval")+
xlab("Continuous HAPS Cutoff")+
geom_vline(xintercept=c(10),linetype=1.2,
color="grey",size=2)
}

plt<-plot_grid(li[[1]],li[[2]],nrow=1,ncol=2)
ggsave("Figure1.D&E.cutoff.pdf",plt,height=5,width=9)


################
#Figure1.D&E down: Survical
custom_theme<-function(){
  theme_survminer()%+replace%
  theme(plot.title=element_text(hjust=0.5))
}

data$group<-rep("",dim(data)[1])
data[which(data$HAPS>10),]$group<-"High HAPS"
data[which(data$HAPS<=10),]$group<-"Low HAPS"


li<-list()

var<-as.data.frame(table(data$Batch))
var
#         Var1 Freq
#1   Training1   30
#2   Training2   34

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

p<-arrange_ggsurvplots(list(li[[1]],li[[2]]),
nrow=1,ncol=2)
ggsave("Figure1.D&E.sur.pdf",p,height = 5,width = 9)

#############################################


all<-read.xlsx("../WES_Cohort.xlsx",sheet=2)
tapply(all$NeoantigenQuality,all$group,summary)

all$NeoantigenQuality<-log10(all$NeoantigenQuality+1)


my_comparisons <- list( c("High HAPS","Low HAPS"))

plt2<-ggplot(all,aes(x=group,y=NeoantigenQuality,fill=group)) +
geom_violin(position = position_dodge(0.9),
alpha = 0.8,width = 1.2,trim = F,color = "grey50") +
geom_boxplot(width = 0.2,
show.legend = F,position = position_dodge(0.9),
color = 'grey20',alpha = 0.8,
outlier.color = 'grey50') +
xlab("")+ylab("Neoantigen Quantily")+
#geom_point(color="black",size=1,position="identity")+
#geom_jitter(position="jitter",colour="grey20")+
theme_classic()+
scale_fill_manual(values = c(c("#CB3425","#3F5688")))+
theme(legend.position="none",
text=element_text(size=15,hjust = 0.5))+
stat_compare_means(comparisons = my_comparisons)

plt1<-ggplot(all,aes(x=NeoantigenQuality,
fill=group))+
geom_density(alpha=0.8)+
theme_classic()+
scale_fill_manual(values = c("#CB3425","#3F5688"))+
xlab("HAPS Group")+
ylab("The Ratio of Neoantigen Quality")+
theme(legend.title=element_blank(),
text=element_text(size=15,hjust = 0.5),
legend.position="top")

best1<-plot_grid(plt1,plt2,nrow=2,ncol=1,rel_heights = c(1,1))

ggsave("Figure1.Quality.pdf",best1,height = 9,width =6)



###################END##########################

