##Figure6
source("/users/zhuxl/HAPS/Code_Script/0.suppressMessages.R")
setwd("/users/zhuxl/HAPS/Code_Script/")

data<-read.xlsx("./Model.xlsx",sheet=1)
p1<-data %>% 
 # mutate_each(funs(rescale), -Sample) %>%
  ggradar(base.size=15,group.point.size=0,
  legend.position="none",
  #plot.title=var[i,1],
  legend.text.size=6)
ggsave("Figure6A.ggradar.pdf",p1,height=4,width=4)


###########
data<-read.xlsx("./Model.xlsx",sheet=2)
data<-data[order(data$Cohort,decreasing=T),]
inTrain = c(1:52)

table(data$Cohort)
table(data$RECIST,data$Response)
data$Cohort
newdata3<-data[,c("HAPS","TMB","Diversity")]
newdata3<-scale(newdata3,scale=TRUE,center=FALSE)
newdata3<-cbind(as.data.frame(data$Response),newdata3)
colnames(newdata3)[1]<-"Response"
newdata3$Response<-as.factor(newdata3$Response)
mdrrClass<-as.factor(newdata3$Response)

unique(data[inTrain,]$Cohort)
trainx = newdata3[inTrain,]
testx = newdata3[-inTrain,]
trainy = mdrrClass[inTrain]
testy = mdrrClass[-inTrain]
dim(trainx);dim(testx)

regr_nn <- train(Response~. , data = trainx, 
method = "nnet")
models<-list(nnet=regr_nn)
predValues = extractPrediction(models,testX = testx[,-1], testY = testy)
head(predValues)
probValues = extractProb(models,testX = testx[,-1], testY = testy)
head(probValues)
table(probValues$dataType)
testProbs = subset(probValues, dataType == "Test")
prob1 = subset(probValues, 
model == "nnet")
library(ROCR)
prob1$label=ifelse(prob1$obs=='DCR',yes=1,0)
pred1 = prediction(prob1$DCR,prob1$label)
perf1 = performance(pred1, measure="tpr", x.measure="fpr")

#prob1<-read.table("clipboard",header=T,sep="\t")
roc.final <- roc(prob1$label, prob1$DCR)
R<-roc.final$auc
CI<-ci.auc(roc.final)
st<-cbind(R,CI)
pdf("Model_ROC.pdf")
plot(roc.final,col="red",
main=paste0("AUC = ",round(as.numeric(CI),3)[2],
"; 95% CI = ",round(as.numeric(CI),3)[1]," - ",
round(as.numeric(CI),3)[3]))
dev.off()

rownames(data)==rownames(prob1)
res<-cbind(data,prob1)
write.csv(res,"final_res.csv",row.names=F)

######################
training<-res[res$Cohort=="Wang-Panel-B",]
cut<-seq(min(training$DCR)+0.001,
max(training$DCR)-0.001,0.001)
H<-data.frame()
for(i in 1:length(cut)){
cutoff<-cut[i]
training$group<-rep("",dim(training)[1])
training[which(training$DCR>cutoff),]$group<-"High Score"
training[which(training$DCR<=cutoff),]$group<-"Low Score"
ma<-as.data.frame(table(training$group,training$Response))
sencitivity<-ma[1,3]/sum(ma[c(1,2),3])
specificity<-ma[4,3]/sum(ma[c(3,4),3])
Youden<-sencitivity+specificity
accuracy<-sum(ma[c(1,4),3])/sum(ma[,3])
a<-cbind(cutoff,sencitivity,specificity,Youden,accuracy)
H<-rbind(H,a)
}
H<-H[order(H$Youden,decreasing=T),]
write.csv(H,"Model_cutoff.csv",row.names=F)
H_cutoff<-H
head(H_cutoff)
dim(H_cutoff)
cutoff<-round(H$cutoff[1],3)
H<-H[which(H$Youden!=2),]
dim(H)
head(H)


p_cutoff<-ggplot(H,aes(x=1-specificity,
y=sencitivity))+
geom_line(na.rm = TRUE,col="#4A4A4A",size=1.5)+
theme_classic()+
ggtitle("")+xlim(c(0,1))+
ylim(c(0,1))+
theme(legend.position="none",
plot.title = element_text(hjust = 0.5),
text=element_text(size=13,hjust = 0.5))+
geom_hline(yintercept=H$sencitivity[1],
linetype=2,color="#F8766D",size=1.5)+
geom_vline(xintercept=1-H$specificity[1],
linetype=2,color="#00BFC4",size=1.5)+
annotate("text", x = 0.5, y = 0.5, 
colour = "red",
label = paste0("Sencitivity = ",round(H$sencitivity[1],3),"\n",
"Specificity = ",round(H$specificity[1],3),"\n",
"Youden = ",round(H$Youden[1],3),"\n",
"Cutoff = ",round(H$cutoff[1],3)))

pdf("Model_cutoff.pdf",height=6,width=6)
p_cutoff
dev.off()

######################
custom_theme<-function(){
  theme_survminer()%+replace%
  theme(plot.title=element_text(hjust=0.5))
}

var<-as.data.frame(table(res$Cohort))

H<-data.frame()
OS<-list()
for(i in 1:dim(var)[1]){
data1<-res[res$Cohort==as.vector(var[i,1]),]
my.surv <- Surv(data1$OS_Months, data1$OS_Event)
fit <- survfit(my.surv ~ data1$group)
data1.survdiff <- survdiff(my.surv ~ data1$group) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))
al<-cbind(p.val,HR,low95,up95)
a<-t(as.data.frame(tapply(data1$OS_Months,data1$group,median)))
b<-cbind(paste0("OS ",as.vector(var[i,1])),a,al)
H<-rbind(H,b)

OS[[i]]<-ggsurvplot(fit, data =data1 ,
           #surv.median.line = "hv",
           #ggtheme = theme(plot.title=element_text(hjust=0.5)),
           ggtheme=custom_theme(),
           conf.int = F,
           conf.int.style = "step",
           censor = T,size=1.2,censor.size=8,
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
           legend.labs=c("High Score","Low Score"),
           font.x = 0,font.y = 13,font.tickslab=13,
           pval = paste(
           paste("P = ",round(p.val,3), sep = ""),
           paste("HR = ",round(HR,2),sep = ""), 
           paste("95%CI: ", paste(round(low95,2), 
           round(up95,2), sep = " - "), sep = ""),
           sep = "\n"),
           pval.coord=c(-2, 0.15))
}

PFS<-list()
for(i in 1:dim(var)[1]){
data1<-res[res$Cohort==as.vector(var[i,1]),]
my.surv <- Surv(data1$PFS_Months, data1$PFS_Event)
fit <- survfit(my.surv ~ data1$group)
data1.survdiff <- survdiff(my.surv ~ data1$group) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))
al<-cbind(p.val,HR,low95,up95)
a<-t(as.data.frame(tapply(data1$PFS_Months,data1$group,median)))
b<-cbind(paste0("PFS ",as.vector(var[i,1])),a,al)
H<-rbind(H,b)

PFS[[i]]<-ggsurvplot(fit, data =data1 ,
           #surv.median.line = "hv",
           #ggtheme = theme(plot.title=element_text(hjust=0.5)),
           ggtheme=custom_theme(),
           conf.int = F,
           conf.int.style = "step",
           censor = T,size=1.2,censor.size=8,
           title=as.vector(var[i,1]),
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
           legend.labs=c("High Score","Low Score"),
           font.x = 0,font.y = 13,font.tickslab=13,
           pval = paste(
           paste("P = ",round(p.val,3), sep = ""),
           paste("HR = ",round(HR,2),sep = ""), 
           paste("95%CI: ", paste(round(low95,2), 
           round(up95,2), sep = " - "), sep = ""),
           sep = "\n"),
           pval.coord=c(-2, 0.15))
}
#H
write.csv(H,"stat_model_time.csv",row.names=F)

p_time<-arrange_ggsurvplots(list(PFS[[1]],PFS[[2]],
OS[[1]],OS[[2]]),
nrow=2,ncol=2)

dev.off()
dev.off()
ggsave("Fig 6 Model_PFS_OS.pdf",p_time,height = 8,width = 10)

######################
data<-read.table("clipboard",header=T,sep="\t")
st<-wilcox.test(DCR~ORR, data=data)
pval<-st$p.value
pval
table(data$RECIST,data$ORR)

st<-wilcox.test(DCR~Response, data=data)
pval<-st$p.value
pval
table(data$RECIST,data$Response)

data<-data[order(data$DCR),]
data<-data[order(data$Cohort),]
data$Score<-data$DCR
data1<-data[data$Cohort=="Wang-Panel-T",]
data1<-data1[order(data1$Score),]
T_label<-c(setdiff(data2$Sample,data1$Sample),data1$Sample)

table(data1$group,data1$ORR)
table(data2$group,data2$ORR)

data2<-data[data$Cohort=="Wang-Panel-B",]
data2$Score<- -(data2$Score)
table(data1$group,data1$ORR)
table(data2$group,data2$ORR)
cutoff<-unique(data$cutoff)
data3<-rbind(data1,data2)


p_fall1<-ggplot(data3,aes(x=Sample,y=Score,fill=RECIST))+
geom_bar(stat="identity",position="dodge")+
theme_bw()+coord_flip()+
xlab("")+ylab("Model Predicition")+
scale_y_continuous(expand=c(0,0))+
scale_x_discrete(limits=as.vector(data2$Sample),
labels=as.vector(data2$Sample))+
theme(legend.position=c(0.15,0.2),
#legend.position="none",
plot.title = element_text(size=15,hjust = 0.5),
axis.text.y=element_blank(),
axis.ticks.y=element_blank())+
guides(fill=guide_legend(reverse=TRUE))+
scale_fill_manual(values = c("#E64B35","#3C5488","#4DBBD5"))+
geom_hline(yintercept=-cutoff,linetype=2,color="grey50")+
geom_hline(yintercept=0,linetype=1,color="black")+
geom_hline(yintercept=cutoff,linetype=2,color="grey50")+
ggtitle("Wang-Panel B cohort & Wang-Panel T cohort")


p_fall2<-ggplot(data3,aes(x=Sample,y=Score,fill=RECIST))+
geom_bar(stat="identity",position="dodge")+
theme_bw()+coord_flip()+
xlab("")+ylab("Model Predicition")+
scale_y_continuous(expand=c(0,0))+
scale_x_discrete(limits=as.vector(T_label),
labels=as.vector(T_label))+
theme(legend.position=c(0.15,0.2),
#legend.position="none",
plot.title = element_text(size=15,hjust = 0.5),
axis.text.y=element_blank(),
axis.ticks.y=element_blank())+
guides(fill=guide_legend(reverse=TRUE))+
scale_fill_manual(values = c("#E64B35","#3C5488","#4DBBD5"))+
geom_hline(yintercept=-cutoff,linetype=2,color="grey50")+
geom_hline(yintercept=0,linetype=1,color="black")+
geom_hline(yintercept=cutoff,linetype=2,color="grey50")+
ggtitle("Wang-Panel B cohort & Wang-Panel T cohort")

p_comb<-plot_grid(p_fall1,p_fall2)

pdf("Model_comb.pdf",height=7,width=14)
p_comb
dev.off()

##########