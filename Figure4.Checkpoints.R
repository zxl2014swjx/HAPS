##Figure4.BCDE
source("/users/zhuxl/HAPS/Code_Script/0.suppressMessages.R")
setwd("/users/zhuxl/HAPS/Code_Script/")

data<-read.xlsx("./Immune_Score.xlsx",sheet=2)
data$group<-rep("",dim(data)[1])
data[which(data$HAPS>10),]$group<-"High HAPS"
data[which(data$HAPS<=10),]$group<-"Low HAPS"
data<-data[,c(1:4,16,5:15)]
dim(data)
#[1] 455  16

data1<-data
m<-as.numeric(data1$HAPS)
n<-as.numeric(data1$pDC)
r<-cor(m,n,method="spearman")
p<-cor.test(m,n)$p.value
r;p

p1<-ggplot(data1,aes(x=HAPS,y=pDC))+
geom_point(na.rm = TRUE,col="#4A4A4A")+
stat_smooth(method =  lm,col="red")+
geom_text(label=NA) + 
annotate("text", x = 20, y = -0.25, size = 4, colour = "red",
label = paste0("Spearman\n","R = ",round(r,3),"\nP = ",
format(p,scientific=TRUE,digit=3)))+
theme(plot.title = element_text(hjust = 0.5))+
theme_bw()

#########################


data1<-melt(data,id.vars=c(1:6))
st<-as.data.frame(table(data1$variable))
li<-list()
my_comparisons<-list(c("High HAPS","Low HAPS"))
for(i in 1:dim(st)[1]){
need<-data1[data1$variable==st[i,1],]

li[[i]]<-ggplot(need,aes(x = group, y = value,
fill = group)) +
geom_violin(#draw_quantiles = c(0.25, 0.5, 0.75),
position = position_dodge(0.9),
alpha = 0.8,width = 1,trim = F,color = "grey50") +
geom_boxplot(width = 0.3,notch=F,
show.legend = F,position = position_dodge(0.9),
color = 'grey20',alpha = 0.8,
outlier.color = 'grey50') +
xlab("")+ylab("immune score")+
facet_wrap(.~variable)+
theme_bw()+
#geom_jitter(position="jitter",colour="grey20")+
scale_fill_manual(values=c("#CB3425","#3F5688"))+
stat_compare_means(comparisons = my_comparisons)+
theme(legend.position="none",
text=element_text(size=13,hjust = 0.5),
plot.title = element_text(size=13,hjust = 0.5))

}

b1<-plot_grid(li[[9]],li[[10]],nrow=1,ncol=2,labels=LETTERS[2])
b2<-plot_grid(p1,nrow=1,ncol=1,labels=LETTERS[3])
b3<-plot_grid(li[[7]],li[[8]],nrow=1,ncol=2,labels=LETTERS[4])
b4<-plot_grid(li[[1]],li[[2]],li[[3]],li[[4]],li[[5]],li[[6]],
nrow=1,ncol=6,labels=LETTERS[5])
best1<-plot_grid(b1,b2,b3,nrow=1,ncol=3)
best2<-plot_grid(best1,b4,nrow=2,ncol=1)

ggsave("Figure 4 down1.pdf",best2,height = 7,width =15)


