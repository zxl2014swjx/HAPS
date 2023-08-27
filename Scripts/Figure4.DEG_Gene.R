library(ggplot2)
library(DESeq2)
library(edgeR)
library(limma)
library(S4Vectors)
library(stats4)
library(BiocGenerics)
library(parallel)
library(BiocGenerics)
library(S4Vectors)
library(IRanges)
library(GenomicRanges)
library(GenomeInfoDb)
library(SummarizedExperiment)
library(Biobase)
library(DelayedArray)
library(matrixStats)
library(BiocParallel)
library(DelayedArray)
library(DESeq2)
library(ggrepel)

exp<-read.table("HAPS.txt",header=T,sep="\t",row.names=1)
dim(exp)
exp[1:4,1:4]
#colnames(exp)[1]<-"Hugo_Symbol"

coding<-read.table("Homo_sapiens.GRCh37.82.gtf.protein_coding.gene.name",header=F,sep="\t")
head(coding)

DEG<-read.csv("DEG_gene.csv")

coding<-as.data.frame(gsub(" ","",coding[,1]))
colnames(coding)<-"Gene"
need<-merge(coding,DEG,by="Gene")
dim(need)

gene<-intersect(rownames(exp),gsub(" ","",coding[,1]))
length(gene)
exp<-exp[gene,]
write.table(exp,"Expression.tsv",sep="\t",quote=F)

group<-read.xlsx("group.xlsx")
length(intersect(colnames(exp),group$Sample))
table(group$Cohort)

Allen<-group[group$Cohort=="Allen",]
exp_Allen<-exp[,Allen$Sample]
scale_exp_Allen<-scale(exp_Allen)

Hugo<-group[group$Cohort=="Hugo",]
exp_Hugo<-exp[,Hugo$Sample]
scale_exp_Hugo<-scale(exp_Hugo)

LUAD<-group[group$Cohort=="LUAD",]
exp_LUAD<-exp[,LUAD$Sample]
scale_exp_LUAD<-scale(exp_LUAD)

LUSC<-group[group$Cohort=="LUSC",]
exp_LUSC<-exp[,LUSC$Sample]
scale_exp_LUSC<-scale(exp_LUSC)

Riaz<-group[group$Cohort=="Riaz",]
exp_Riaz<-exp[,Riaz$Sample]
scale_exp_Riaz<-scale(exp_Riaz)

Snyder<-group[group$Cohort=="Snyder",]
exp_Snyder<-exp[,Snyder$Sample]
scale_exp_Snyder<-scale(exp_Snyder)

#group<-group[group$Cancer=="SKCM",]

High<-as.data.frame(group[group$group=="High HAPS",]$Sample)
colnames(High)<-"Sample"
Low<-as.data.frame(group[group$group=="Low HAPS",]$Sample)
colnames(Low)<-"Sample"
dim(High)
dim(Low)

table(group$Cohort)
scale_exp1<-cbind(scale_exp_Allen,scale_exp_Hugo,
scale_exp_LUAD,scale_exp_LUSC,
scale_exp_Riaz,scale_exp_Snyder)

scale_exp<-scale(exp)
scale_exp[1:4,1:4]
limm(scale_exp1,Low,High,"Low_VS_High",0.05,1)

exp<-scale_exp1
control<-Low
case<-High
sample<-"Low_VS_High"
fdr<-0.05
logfc<-1


limm<-function(exp,control,case,sample,fdr,logfc){
  group1 <- rep('control',dim(control)[1])
  group2 <- rep('case',dim(case)[1])
  grouplist<-c(group1,group2)
  table(grouplist)
  design <- model.matrix(~0+factor(grouplist))
  colnames(design)=levels(factor(grouplist))
  names<-vector()
  samples<-as.vector(rbind(control,case)[,1])
  #for(j in 1:length(samples)){
  #a<-grep(samples[j],colnames(exp))
  #names<-c(names,a)}
  exprSet<-exp[,samples]
  dim(exprSet)
  rownames(design)=colnames(exprSet)

  fit <- lmFit(exprSet, design)
  cont.matrix <- makeContrasts("control-case", levels=design)
  fit2 <- contrasts.fit(fit,cont.matrix)
  fit2 <- eBayes(fit2)
  tempOutput<-topTable(fit2,coef=1, n=Inf)

  DEG_all<-data.frame(as.vector(rownames(tempOutput)),tempOutput$logFC,tempOutput$AveExpr,tempOutput$P.Value,tempOutput$adj.P.Val,stringsAsFactors = FALSE)
  colnames(DEG_all)<-c("Gene","logFC","AveExpr","pval","adjp")

  up_index<-which(DEG_all$pval <= fdr & DEG_all$logFC >= logfc)
  DEG_up<-DEG_all[up_index,]
  DEG_up$Deg<-"up"
  down_index<-which(DEG_all$pval <= fdr & DEG_all$logFC <= (-logfc))
  DEG_down<-DEG_all[down_index,]
  DEG_down$Deg<-"down"
  DEG_select<-rbind(DEG_up,DEG_down)
  dim(DEG_select)

  name1<-paste0(sample,".txt")
  name2<-paste0(sample,"_select.txt")
  write.table(DEG_all,name1,sep="\t",quote=FALSE,row.names =F)
  write.table(DEG_select,name2,sep="\t",quote=FALSE,row.names =F)
}


#######################


data<-read.table("Low_VS_High_all.txt",header=T,sep="\t")
data<-read.xlsx("DEG_gene.xlsx")
head(data)
data$log10FC<- log10(2^data$logFC+1)


data$FOLD<-rep("unchange",dim(data)[1])
data[which(data$logFC>1 & data$pval < 0.05),]$FOLD<- "up"
data[which(data$logFC< -1 & data$pval < 0.05),]$FOLD<- "down"
table(data$FOLD)
write.xlsx(data,"DEG_gene1.xlsx")

ta<-as.data.frame(table(data$FOLD))
up_num<-ta[which(ta[,1]=="up"),2]
down_num<-ta[which(ta[,1]=="down"),2]
unchange_num<-ta[which(ta[,1]=="unchange"),2]
data$FOLD<-as.factor(data$FOLD)
table(data$FOLD)

p1<-ggplot(data,mapping=aes(x=logFC, y=-1*log10(adjp), color=FOLD))+
 geom_point(size=1)+
 xlab("log2(FC)")+ylab("-log10(FDR)")+
 ggtitle("Low VS High HAPS DEG")+
 theme(plot.title = element_text(hjust = 0.5))+
 theme(plot.title=element_text(face="bold",size=14))+
 scale_color_manual(values=c("green","grey","red"),
 labels=c(paste0("down: ",down_num),paste0("unchange: ",unchange_num),paste0("up: ",up_num)))

ggsave("ALL_DEG.pdf",p1,height=5,width=7)

########################

data_posi<-data[which(data$logFC>1 & data$pval < 0.05),]
data_nege<-data[which(data$logFC< -1 & data$pval < 0.05),]
data_un<-data[-which((data$logFC< -1 | data$logFC>1)& data$pval < 0.05),]
dim(data_posi);dim(data_nege);dim(data_un)

p<-ggplot()+
geom_point(data,mapping=aes(x=logFC, y=-log10(pval),
color=FOLD))+
geom_vline(xintercept=c(-1,0,1),color="grey",linetype=2)+
geom_hline(yintercept=-log10(0.05),color="grey",linetype=2)+
xlab("log2(FC)")+ylab("-log10(FDR)")+
theme_classic()+
theme(panel.border=element_blank(),
legend.position="top",
axis.line.x=element_line(color="black"))+
geom_point(aes(x=data_posi$logFC,
y=-log10(data_posi$pval)),color="pink")+
geom_point(aes(x=data_nege$logFC,
y=-log10(data_nege$pval)),color="lightblue")+
geom_point(aes(x=data_un$logFC,
y=-log10(data_un$pval)),color="grey")+

scale_color_manual(values=c("lightblue","grey","pink"),
labels=c(paste0("down: ",down_num),paste0("unchange: ",unchange_num),paste0("up: ",up_num)))+

geom_label_repel(data=data_posi,aes(x=logFC, 
y=-log10(data_posi$pval),
label =as.vector(data_posi$Gene)),
point.padding=unit(0.5, "lines"),
segment.colour = "grey",fill="pink", size=3,
fontface="bold",,color="black",
arrow = arrow(length=unit(0.01,"npc")))+

geom_label_repel(data=data_nege,aes(x=logFC, 
y=-log10(data_nege$pval),
label =as.vector(data_nege$Gene)),
point.padding=unit(0.5, "lines"),
segment.colour = "grey",fill="lightblue",size=3,
fontface="bold",color="black",
arrow = arrow(length=unit(0.01,"npc")))



pdf("Volcano_0.05_1.pdf",height=5,width=8)
p
dev.off()


