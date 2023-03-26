##Figure5.A
source("/users/zhuxl/HAPS/Code_Script/0.suppressMessages.R")
setwd("/users/zhuxl/HAPS/Code_Script/Figure 5")

data<-read.xlsx("./TCR.xlsx",sheet=1)

data1<-data
m<-as.numeric(data1$HAPS)
n<-as.numeric(data1$Diversity)
r<-cor(m,n,method="spearman")
p<-cor.test(m,n)$p.value
r;p

p1<-ggplot(data1,aes(x=HAPS,y=Diversity))+
geom_point(na.rm = TRUE,col="#4A4A4A")+
stat_smooth(method =  lm,col="red")+
geom_text(label=NA) + 
annotate("text", x = 5, y = 5, size = 4, colour = "red",
label = paste0("Spearman\n","R = ",round(r,3),"\nP = ",
format(p,scientific=TRUE,digit=3)))+
theme(plot.title = element_text(hjust = 0.5))+
theme_bw()

ggsave("Figure5.HAPS_Diversity.pdf",p1,height=4,width=4)

###############
data<-read.xlsx("./Wang_Panel.xlsx",sheet=1)
data1<-data[data$Cohort=="Wang-Panel-T",]
my.surv <- Surv(data1$OS_Months, data1$OS_Event)
fit <- survfit(my.surv ~ data1$HAPS_Diversity)
data1.survdiff <- survdiff(my.surv ~ data1$HAPS_Diversity) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))

T_D<-ggsurvplot(fit, data =data1 ,
           ggtheme=custom_theme(),
           conf.int = F,
           conf.int.style = "step",
           censor = T,size=2,censor.size=8,
           title="Wang-Panel-T HAPS & Diversity",
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
           legend.labs=c("High HAPS & High Diversity",
           "High HAPS & Low Diversity",
           "Low HAPS & High Diversity",
           "Low HAPS & Low Diversity"),
           font.x = 0,font.y = 13,font.tickslab=13,
           pval = paste(
           paste("P = ",round(p.val,4), sep = ""),
           paste("HR = ",round(HR,2),sep = ""), 
           paste("95%CI: ", paste(round(low95,2), 
           round(up95,2), sep = " - "), sep = ""),
           sep = "\n"),
           pval.coord=c(-2, 0.15))

###############

data<-read.xlsx("./Wang_Panel.xlsx",sheet=1)
data1<-data[data$Cohort=="Wang-Panel-B",]

my.surv <- Surv(data1$OS_Months, data1$OS_Event)
fit <- survfit(my.surv ~ data1$HAPS_Diversity)
data1.survdiff <- survdiff(my.surv ~ data1$HAPS_Diversity) 
p.val = 1 - pchisq(data1.survdiff$chisq, length(data1.survdiff$n) - 1)
p.val
HR = (data1.survdiff$obs[1]/data1.survdiff$exp[1])/(data1.survdiff$obs[2]/data1.survdiff$exp[2])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1])) 
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data1.survdiff$exp[2]+1/data1.survdiff$exp[1]))

B_D<-ggsurvplot(fit, data =data1 ,
           #surv.median.line = "hv",
           #ggtheme = theme(plot.title=element_text(hjust=0.5)),
           ggtheme=custom_theme(),
           conf.int = F,
           conf.int.style = "step",
           censor = T,size=2,censor.size=8,
           title="Wang-Panel-B HAPS & Diversity",
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
           legend.labs=c("High HAPS & High Diversity",
           "High HAPS & Low Diversity",
           "Low HAPS & High Diversity",
           "Low HAPS & Low Diversity"),
           font.x = 0,font.y = 13,font.tickslab=13,
           pval = paste(
           paste("P = ",round(p.val,4), sep = ""),
           paste("HR = ",round(HR,2),sep = ""), 
           paste("95%CI: ", paste(round(low95,2), 
           round(up95,2), sep = " - "), sep = ""),
           sep = "\n"),
           pval.coord=c(-2, 0.15))


p<-arrange_ggsurvplots(list(T_D,B_D),nrow=1,ncol=2)
ggsave("Figure5.Survival.pdf",p,height = 5,width = 12)


###############


data<-read.xlsx("./TCR.xlsx",sheet=1)


data$Delta<-data$Post_Clonality-data$Baseline_Clonality
st<-compare_means(Delta~benefit, data=data)
p=st$p.adj

data1<-melt(data,id.vars=c(1:4))
p1<-ggplot(data1,aes(x=variable,y=value,
col=clonality_cut,group=Pat))+
geom_point()+
ggtitle(paste0("Wilcox Test P =",p))+
geom_line()+xlab("Treatment")+ylab("Clonality")+
facet_wrap(.~benefit,nrow=1,ncol=2)+
theme_bw()+
theme(legend.position="none",
text=element_text(size=12,hjust = 0.5),
plot.title = element_text(size=12,hjust = 0.5))+
scale_color_manual(values = c("#3F5688","#CB3425"))

ggsave("Figure5D.ggpaired.pdf",p1,height=5,width=5)


#############
data<-read.xlsx("TCR.xlsx",sheet=2)

my_comparisons1 <- list(c("High HAPS","Low HAPS"))

my_comparisons2 <- list(
c("High HAPS & High Diversity","High HAPS & Low Diversity"),
c("High HAPS & High Diversity","Low HAPS & High Diversity"),
c("High HAPS & High Diversity","Low HAPS & Low Diversity"),
c("High HAPS & Low Diversity","Low HAPS & High Diversity"),
c("High HAPS & Low Diversity","Low HAPS & Low Diversity"),
c("Low HAPS & High Diversity","Low HAPS & Low Diversity")
)


my_comparisons3 <- list(
c("High HAPS & Post_tre","High HAPS & Baseline"),
c("High HAPS & Post_tre","Low HAPS & Post_tre"),
c("High HAPS & Post_tre","Low HAPS & Baseline"),
c("High HAPS & Baseline","Low HAPS & Post_tre"),
c("High HAPS & Baseline","Low HAPS & Baseline"),
c("Low HAPS & Post_tre","Low HAPS & Baseline")
)


plt1<-ggviolin(data, x="Group", 
y="edit_distance",#bxp.errorbar=TRUE,
fill = "Group", 
add=c("mean","boxplot","points"),
xlab="",ylab="Edit Distance",
palette="npg",
trim=T,
order=c("High HAPS","Low HAPS")
)+
#coord_flip()+
theme_classic()+
stat_compare_means(comparisons = my_comparisons1)+
#stat_compare_means(label = "p.signif", method = "",ref.group = ".all.")+
#stat_compare_means(label.y =20)+
theme(legend.position='none',
axis.text.x = element_text(angle = 45,
hjust = 1,color = 'black'),
plot.title = element_text(size=15,hjust = 0.5))+
scale_y_break(c(6, 7),scales="free")+
scale_y_break(c(11,12),scales="free")


plt2<-ggviolin(data, x="label1", 
y="edit_distance",#bxp.errorbar=TRUE,
fill = "label1", 
add=c("mean","boxplot","points"),
xlab="",ylab="Edit Distance",
palette="npg",
trim=T,
order=c("High HAPS & High Diversity",
"High HAPS & Low Diversity",
"Low HAPS & High Diversity",
"Low HAPS & Low Diversity")
)+
#coord_flip()+
theme_classic()+
stat_compare_means(comparisons = my_comparisons2)+
#stat_compare_means(label = "p.signif", method = "",ref.group = ".all.")+
#stat_compare_means(label.y =20)+
theme(legend.position='none',
axis.text.x = element_text(angle = 45,
hjust = 1,color = 'black'),
plot.title = element_text(size=15,hjust = 0.5))+
scale_y_break(c(6, 7),scales="free")+
scale_y_break(c(11,12),scales="free")



plt3<-ggviolin(data, x="label2", 
y="edit_distance",#bxp.errorbar=TRUE,
fill = "label2", 
add=c("mean","boxplot","points"),
xlab="",ylab="Edit Distance",
palette="npg",
trim=T,
order=c("High HAPS & Baseline","High HAPS & Post_tre",
"Low HAPS & Baseline","Low HAPS & Post_tre")
)+
#coord_flip()+
theme_classic()+
stat_compare_means(comparisons = my_comparisons3)+
#stat_compare_means(label = "p.signif", method = "",ref.group = ".all.")+
#stat_compare_means(label.y =20)+
theme(legend.position='none',
axis.text.x = element_text(angle = 45,
hjust = 1,color = 'black'),
plot.title = element_text(size=15,hjust = 0.5))+
scale_y_break(c(6, 7),scales="free")+
scale_y_break(c(11,12),scales="free")

al<-plot_grid(plt1,plt2,plt3,nrow=1,ncol=3,
rel_widths = c(1,2,2))

ggsave("Figure5.Edit Distance.pdf",al,height=5,width=13)
#######################


