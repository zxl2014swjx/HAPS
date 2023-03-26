##Figure4.FG
source("/users/zhuxl/HAPS/Code_Script/0.suppressMessages.R")
setwd("/users/zhuxl/HAPS/Code_Script/Figure 4")

data<-read.xlsx("./Immune_Score.xlsx",sheet=3)

data_posi<-data[which(data$log10FC>1 & data$adjp < 0.05),]
data_nege<-data[which(data$log10FC< -1 & data$adjp < 0.05),]
data_un<-data[-which((data$log10FC< -1 | data$log10FC>1)& data$adjp < 0.05),]
dim(data_posi);dim(data_nege);dim(data_un)
#[1] 15  7
#[1] 22  7
#[1] 45  7

p<-ggplot()+
geom_point(data,mapping=aes(x=log10FC, y=-log10(adjp),
color=Group))+
geom_vline(xintercept=c(-1,0,1),color="grey",linetype=2)+
geom_hline(yintercept=-log10(0.05),color="grey",linetype=2)+
xlab("log10(FC)")+ylab("-log10(FDR)")+
theme_classic()+
theme(panel.border=element_blank(),
legend.position="top",
axis.line.x=element_line(color="black"))+
geom_point(aes(x=data_posi$log10FC,
y=-log10(data_posi$adjp)),color="pink")+
geom_point(aes(x=data_nege$log10FC,
y=-log10(data_nege$adjp)),color="lightblue")+
geom_point(aes(x=data_un$log10FC,
y=-log10(data_un$adjp)),color="grey")+
scale_color_manual(values=c("lightblue","grey","pink"),
labels=c(paste0("Low HAPS"),
paste0("unchange"),
paste0("High HAPS")))+
geom_label_repel(data=data_posi,aes(x=log10FC, 
y=-log10(data_posi$adjp),
label =as.vector(data_posi$Gene)),
point.padding=unit(0.5, "lines"),
segment.colour = "grey",fill="pink", size=3,
fontface="bold",,color="black",
arrow = arrow(length=unit(0.01,"npc")))+
geom_label_repel(data=data_nege,aes(x=log10FC, 
y=-log10(data_nege$adjp),
label =as.vector(data_nege$Gene)),
point.padding=unit(0.5, "lines"),
segment.colour = "grey",fill="lightblue",size=3,
fontface="bold",color="black",
arrow = arrow(length=unit(0.01,"npc")))

ggsave("Volcano1.pdf",p,height=6,width=6)


##################################
data<-read.xlsx("./Immune_Score.xlsx",sheet=4)
data<-data[order(data$pvalue,decreasing=TRUE),]
data<-data[order(data$ID,decreasing=TRUE),]

p<-ggplot(data,aes(x=log10(1/pvalue+1),
y=Description))+
geom_point(aes(size=rank,col=log10(1/pvalue+1)))+
theme_bw()+
scale_colour_gradient(low = "blue",high = "red",na.value = "white")+
scale_y_discrete(limits=as.vector(data$Description),labels=as.vector(data$Description))+
geom_vline(xintercept=log10(1/0.5+1),linetype=2,color="grey50")+
geom_hline(yintercept=7.5,linetype=2,color="grey50")+
geom_hline(yintercept=12.5,linetype=2,color="grey50")+
ylab("")+xlab("")+
theme(legend.position="top",
axis.text.y=element_blank(),
axis.ticks.y=element_blank(),
#axis.text.y=element_text(hjust=1,size=5),
text=element_text(hjust = 0.5))+
geom_text(aes(label = Description),col="grey50",
vjust = "inward")

ggsave("Fig 4 sig_pathway.pdf",p,height=6,width=10)

#################END####################
