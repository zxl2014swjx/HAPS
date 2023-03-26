##Figure3.Sur
source("/users/zhuxl/HAPS/Code_Script/0.suppressMessages.R")
setwd("/users/zhuxl/HAPS/Code_Script/Figure 3")

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

data$group<-rep("",dim(data)[1])
data[which(data$HAPS>10),]$group<-"High HAPS"
data[which(data$HAPS<=10),]$group<-"Low HAPS"

##############
###Figure3A.HAPS Group & LOH

st<-table(data$group,data$LOH)
st1<-as.data.frame(st)
colnames(st1)<-c("group","LOH","Num")
p<-fisher.test(st)$p.value
pval<-format(p,scientific=TRUE,digit=3)
pval<-round(p,3)
p1<-ggplot(st1,aes(x=group,y=Num,fill=LOH))+
geom_bar(stat="identity",position="fill")+
ggtitle(paste0("Fisher's Exact Test\n", 
"p = ",pval))+ylab("Percentage")+
theme_classic()+
scale_fill_manual(values = c("#CB3425","#3F5688"))+
scale_y_continuous(expand=c(0,0))+
theme(#legend.position="top",
plot.title = element_text(size=13,hjust = 0.5))
ggsave("Fig3 Group LOH.pdf",p1,height=4,width=4)

############################
data1<-data
my.surv <- Surv(data1$OS_Months, data1$OS_Event)
fit <- survfit(my.surv ~ data1$group)
data1.survdiff <- survdiff(my.surv ~ data1$group) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))
al<-cbind("ALL",dim(data1)[1],p.val,HR,low95,up95)

data1<-data[which(data$LOH=="LOH=0"),]
my.surv <- Surv(data1$OS_Months, data1$OS_Event)
fit <- survfit(my.surv ~ data1$group)
data1.survdiff <- survdiff(my.surv ~ data1$group) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))
al0<-cbind(paste0("ALL"," LOH=0"),
dim(data1)[1],p.val,HR,low95,up95)

data1<-data[which(data$LOH=="LOH=1"),]
my.surv <- Surv(data1$OS_Months, data1$OS_Event)
fit <- survfit(my.surv ~ data1$group)
data1.survdiff <- survdiff(my.surv ~ data1$group) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))
al1<-cbind(paste0("ALL"," LOH=1"),
dim(data1)[1],p.val,HR,low95,up95)

H0<-rbind(al,al0,al1)
colnames(H0)[1:2]<-c("batch","V2")
###############
H1<-data.frame()
var<-as.data.frame(table(data$Batch))
for(i in 1:dim(var)[1]){
batch<-as.vector(var[i,1])

data1<-data[data$Batch==batch,]
my.surv <- Surv(data1$OS_Months, data1$OS_Event)
fit <- survfit(my.surv ~ data1$group)
data1.survdiff <- survdiff(my.surv ~ data1$group) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))
al<-cbind(batch,dim(data1)[1],p.val,HR,low95,up95)

data1<-data[which(data$Batch==batch &
data$LOH=="LOH=0"),]
my.surv <- Surv(data1$OS_Months, data1$OS_Event)
fit <- survfit(my.surv ~ data1$group)
data1.survdiff <- survdiff(my.surv ~ data1$group) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))
al0<-cbind(paste0(batch," LOH=0"),
dim(data1)[1],p.val,HR,low95,up95)

data1<-data[which(data$Batch==batch &
data$LOH=="LOH=1"),]
my.surv <- Surv(data1$OS_Months, data1$OS_Event)
fit <- survfit(my.surv ~ data1$group)
data1.survdiff <- survdiff(my.surv ~ data1$group) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))
al1<-cbind(paste0(batch," LOH=1"),
dim(data1)[1],p.val,HR,low95,up95)

al2<-rbind(al,al0,al1)
H1<-rbind(H1,al2)

}

H<-rbind(H0,H1)
colnames(H)[1:2]<-c("Batch","SampleSize")
write.csv(H,"Fig3 Cox LOH.csv",row.names=F)

###########################

OS_data<-read.table("clipboard",header=FALSE,sep="\t")
OS_draw<-read.table("clipboard",header=TRUE,sep="\t")

pdf("Figure3.Batch_LOH.pdf",height=6,width=10,onefile=FALSE)
forestplot::forestplot(OS_data, OS_draw,
vertices = TRUE, 
graph.pos = 6,
new_page = FALSE,
hrzl_lines = list(
"2" = gpar(lty=1,columns=c(1:5),col = "#000044"),
"5" = gpar(lty=1,columns=c(1:5),col = "#000044"),
"8" = gpar(lty=1,columns=c(1:5),col = "#000044"),
"11" = gpar(lty=1,columns=c(1:5),col = "#000044"),
"14" = gpar(lty=1,columns=c(1:5),col = "#000044"),
"17" = gpar(lty=1,columns=c(1:5),col = "#000044")
), 
clip=c(0.001,17),xlog=TRUE,
col = fpColors(
    box = c('blue',"black",'blue',"black",'blue',"black",
'blue',"black",'blue',"yellow",'red',
'blue',"black",'blue',"yellow",'red',
'blue',"black",'blue',"yellow",'red'),
    line = "black"))
dev.off()


