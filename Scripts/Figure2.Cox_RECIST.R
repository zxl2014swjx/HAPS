##Figure2.Cox_Stage_RECIST_density
source("/users/zhuxl/HAPS/Code_Script/0.suppressMessages.R")
setwd("/users/zhuxl/HAPS/Code_Script")

data<-read.xlsx("./WES_Cohort.xlsx",sheet=3)
OS <- Surv(data$OS_Months,data$OS_Event==1)

sing_OS<-data.frame()
for(i in 5:10){
a<-coxph(OS~data[,i],data=data)
b<-as.data.frame(summary(a)$coefficients)
c<-as.data.frame(summary(a)$conf.int)
d<-cbind(b,c)
sing_OS<-rbind(sing_OS,d)
}
write.csv(sing_OS,"sing_OS.csv")

mutil_OS<-data.frame()
a<-coxph(OS~TMB+PD_L1+Stage+Age+Gender+group,
data=data)
b<-as.data.frame(summary(a)$coefficients)
c<-as.data.frame(summary(a)$conf.int)
d<-cbind(b,c)
mutil_OS<-rbind(mutil_OS,d)
write.csv(mutil_OS,"mutil_OS.csv")

OS_data<-read.table("clipboard",header=FALSE,sep="\t")
###Multi Cox Analysis	Pvalue	HR	95%CI
OS_draw<-read.table("clipboard",header=TRUE,sep="\t")
###HR	LowerCI	UpperCI


pdf("Fig 2 single_cox.pdf",height=4,width=8)
forestplot::forestplot(OS_data, 
OS_draw,
vertices = TRUE, 
graph.pos = 5,
new_page = FALSE,
hrzl_lines = list("2" = gpar(lty=1,columns=c(1:4),col = "#000044"),
"7" = gpar(lty=1,columns=c(1:4),col = "#000044")), 
clip=c(0.001,17),xlog=TRUE,
col=fpColors(box="blue",line="black",hrz_lines = "red"))
dev.off()

pdf("Fig 2 multi_cox.pdf",height=4,width=8)
forestplot::forestplot(OS_data, 
OS_draw,
vertices = TRUE, 
graph.pos = 5,
new_page = FALSE,
hrzl_lines = list("2" = gpar(lty=1,columns=c(1:4),col = "#000044"),
"7" = gpar(lty=1,columns=c(1:4),col = "#000044")), 
clip=c(0.001,17),xlog=TRUE,
col=fpColors(box="blue",line="black",hrz_lines = "red"))
dev.off()


###############

data<-read.xlsx("./WES_Cohort.xlsx",sheet=1)


my_comparisons<-list(c("DCB","NDB"))

p1<-ggplot(data,aes(x = Response, y = HAPS,
fill = Response)) +
geom_violin(#draw_quantiles = c(0.25, 0.5, 0.75),
position = position_dodge(0.9),
alpha = 0.8,width = 1,trim = F,color = "grey50") +
geom_boxplot(width = 0.3,notch=T,
show.legend = F,position = position_dodge(0.9),
color = 'grey20',alpha = 0.8,
outlier.color = 'grey50') +
xlab("Response")+ylab("HAPS")+
#ggtitle(as.vector(var[i,1]))+
#facet_wrap(.~Cohort)+
theme_bw()+
#geom_jitter(position="jitter",colour="grey20")+
scale_fill_manual(values=c("#CB3425","#3F5688"))+
stat_compare_means(comparisons = my_comparisons)+
theme(legend.position="none",
text=element_text(size=13,hjust = 0.5),
plot.title = element_text(size=13,hjust = 0.5))
ggsave("Figure2.Response_HAPS.pdf",p1,height=5,width=5)

########

data<-read.xlsx("./WES_Cohort.xlsx",sheet=1)
data<-data[which(data$RECIST=="CR" |
data$RECIST=="PD" |
data$RECIST=="PR" |
data$RECIST=="SD"),]
table(data$RECIST,data$group)

data$order<-data$RECIST
data[data$RECIST=="CR",]$order<-"1_CR"
data[data$RECIST=="PR",]$order<-"2_PR"
data[data$RECIST=="SD",]$order<-"3_SD"
data[data$RECIST=="PD",]$order<-"4_PD"
table(data$order,data$group)

data<-data[order(data$HAPS),]
data<-data[order(data$order),]

my_comparisons <- list( 
c("1_CR", "2_PR"),c("2_PR", "3_SD"),c("3_SD", "4_PD"),
c("1_CR", "3_SD"),c("2_PR", "4_PD"),c("1_CR", "4_PD"))

compare_means(HAPS~order, data=data)

p1<-ggplot(data,aes(x=Sample,y=HAPS,fill=order))+
geom_bar(stat="identity",position="dodge")+
theme_classic()+xlab("")+
scale_y_continuous(expand=c(0,0))+
scale_x_discrete(limits=as.vector(data$Sample),
labels=as.vector(data$Sample))+
theme(plot.title = element_text(size=15,hjust = 0.5),
axis.text.x=element_blank(),
axis.ticks.x=element_blank())+
scale_fill_manual(values = c("#E64B35","#4DBBD5",
"#00A087","#3C5488"))+
geom_hline(yintercept=c(10),linetype=2,color="black")

ggsave("Figure2.Sample_RECIST.pdf",p1,height=5,width=5)


###########
data<-read.xlsx("./Wang_Panel.xlsx",sheet=1)

p1<-ggplot(data,aes(x=HAPS,fill=Batch))+
geom_density(alpha=0.8)+
theme_classic()+
scale_fill_manual(values = c("#E64B35","#4DBBD5"))+
xlab("HLA-I Antigen Presentation Score")+
ylab("The Ratio of HAPS")+
theme(legend.title=element_blank(),
text=element_text(size=15,hjust = 0.5),
legend.position="top")
ggsave("Figure2.panel_HAPS_density.pdf",p1,height=6,width=6)




