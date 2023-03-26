##Figure1.A
source("/users/zhuxl/HAPS/Code_Script/0.suppressMessages.R")
setwd("/users/zhuxl/HAPS/Code_Script/Figure 4")

data<-read.xlsx("./Immune_Score.xlsx",sheet=1)
data$group<-rep("",dim(data)[1])
data[which(data$HAPS>10),]$group<-"High HAPS"
data[which(data$HAPS<=10),]$group<-"Low HAPS"
data<-data[,c(1:4,27,5:26)]
dim(data)
#[1] 455  27

allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB",
            "#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3",
            "#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#40E0D0","#5F9EA0","#FF1493",
            "#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
            "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347",
            "#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

n <- 22
g <- erdos.renyi.game(n, 0.5)
la <- layout.circle(g)
radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}
lab.locs <- radian.rescale(x=1:n, direction=-1, start=0)

######################

data1<-data[data$group=="High HAPS",]
data1[1:6,1:6]
data2<-melt(data1,id.vars=c(1:5))
colnames(data2)[6:7]<-c("ImmuneCells","ImmuneScore")
st1<-as.data.frame(tapply(data2$ImmuneScore,
data2$ImmuneCells,sum))
st2<-cbind(rownames(st1),st1)
colnames(st2)<-colnames(data2)[6:7]
rownames(st2)<-1:dim(st2)[1]
write.csv(st2,"High.st.csv")
H<-data.frame()
for(i in 6:dim(data1)[2]){
for(j in 6:dim(data1)[2]){
m<-as.numeric(data1[,i])
n<-as.numeric(data1[,j])
r<-cor(m,n,method="spearman")
p<-cor.test(m,n)$p.value
a<-cbind(colnames(data1)[i],colnames(data1)[j],r,p)
H<-rbind(H,a)
}
}
H$PN<-"Neg"
H[(H$r>0),5]<-"Pos"
H$logp<-  log10(1/as.numeric(as.vector(H$p)))
H$R<-abs(as.numeric(as.vector(H$r)))*10
H<-H[-which(H[,1]==H[,2]),]
H<-H[H$p<0.05,]
write.csv(H,"High.csv",row.names=F)

mynet<-H#[,c(1,2,5,6)]
net <- graph_from_data_frame(mynet, 
directed = TRUE, vertices=NULL) 
karate_groups <- cluster_optimal(net)
coords <- layout_in_circle(net)  
E(net)$width  <- E(net)$logp
E(net)$width  <- E(net)$R
E(net)$color  <- E(net)$PN
E(net)$color[E(net)$PN=="Pos"]<-"pink"
E(net)$color[E(net)$PN=="Neg"]<-"lightblue"

 
pdf("High.pdf",height=8,width=8)
plot(net, 
     edge.arrow.size=0.1, 
     edge.curved=0.1, 
     vertex.label.cex=1,   
     vertex.color=allcolour,
     vertex.label.dist=1,
     vertex.label.degree=lab.locs,
     vertex.size=as.vector(st2[,2]),
     vertex.frame.color=allcolour,
     vertex.label.color="black",
     layout = coords)
dev.off()

###################

data1<-data[data$group=="Low HAPS",]
data1[1:6,1:6]
data2<-melt(data1,id.vars=c(1:5))
colnames(data2)[6:7]<-c("ImmuneCells","ImmuneScore")
st1<-as.data.frame(tapply(data2$ImmuneScore,
data2$ImmuneCells,sum))
st2<-cbind(rownames(st1),st1)
colnames(st2)<-colnames(data2)[6:7]
rownames(st2)<-1:dim(st2)[1]
write.csv(st2,"Low.st.csv")

H<-data.frame()
for(i in 6:dim(data1)[2]){
for(j in 6:dim(data1)[2]){
m<-as.numeric(data1[,i])
n<-as.numeric(data1[,j])
r<-cor(m,n,method="spearman")
p<-cor.test(m,n)$p.value
a<-cbind(colnames(data1)[i],colnames(data1)[j],r,p)
H<-rbind(H,a)
}
}
H$PN<-"Neg"
H[(H$r>0),5]<-"Pos"
H$logp<-  log10(1/as.numeric(as.vector(H$p)))
H$R<-abs(as.numeric(as.vector(H$r)))*10
H<-H[-which(H[,1]==H[,2]),]
H<-H[H$p<0.05,]
write.csv(H,"Low.csv",row.names=F)

mynet<-H#[,c(1,2,5,6)]
net <- graph_from_data_frame(mynet, 
directed = TRUE, vertices=NULL) 
karate_groups <- cluster_optimal(net)
coords <- layout_in_circle(net)  
E(net)$width  <- E(net)$logp
E(net)$width  <- E(net)$R
E(net)$color  <- E(net)$PN
E(net)$color[E(net)$PN=="Pos"]<-"pink"
E(net)$color[E(net)$PN=="Neg"]<-"lightblue"

pdf("Low.pdf",height=8,width=8)
plot(net, 
     edge.arrow.size=0.1, 
     edge.curved=0.1, 
     vertex.label.cex=1,   
     vertex.color=allcolour,
     vertex.label.dist=1,
     vertex.label.degree=lab.locs,
     vertex.size=as.vector(st2[,2]),
     vertex.frame.color=allcolour,
     vertex.label.color="black",
     layout = coords)
dev.off()





