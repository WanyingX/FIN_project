#!/usr/bin/env Rscript

# This script is used to calculate loop size under each condition.

args = commandArgs(trailingOnly=TRUE)
k=as.numeric(args[1])
cut=as.numeric(args[2])
filename=args[3]
prefix=args[4]

#cutoff=as.numeric(args[4])
library(ggplot2)
library(MASS)
#library(viridis)
df=read.table(filename,sep='\t',stringsAsFactors=F,header=T,row.names=1)
df$fc=df$KO/df$WT
a=df[order(-df$WT),]
xintercept=min(a[1:500000,1])
yintercept=xintercept*k
cutoff=k*xintercept


gain=df[df$fc>4,]
lost=df[df$fc<.25,]
retain=df[which(abs(df$residuals)<cutoff&df$WT>xintercept&df$KO>yintercept),]
#lost=df[which(df$residuals<=(-cutoff)),]
#loop=list(WT=wt,Loss=lost,KO=ko,Gain=gain)
loop=list(lost=lost,retain=retain,gain=gain)
#Type=c('WT','Loss','KO','Gain')
Type=names(loop)
xlabel=c()
for(i in Type){
	dist=loop[[i]]$dist
	dist=dist[which(dist<(cut*1000))]
	n=formatC(median(dist),format="f", big.mark = ",", digits=0)
#        n=formatC(dim(loop[[i]])[1],format="f", big.mark = ",", digits=0)
        xlabel=c(xlabel,paste(i,n,sep='\n'))
}
names(xlabel)=c("Lost","Retain","Gain")
lost$Type="Lost"
gain$Type="Gain"
retain$Type="Retain"
df1=rbind(gain,lost,retain)
#col=c("white", "blue",colorRampPalette(c("orange","red"))(8))
if(max(df1$dist>30000)){df1$dist=df1$dist/1000}
#col=c("#999999","blue","orange","red")
df1$Type=factor(df1$Type,levels=c("Lost","Retain","Gain"))
data=df1
data=data[which(data$dist<cut),]
write.table(df1,'dist.bed',sep='\t',quote=F,row.names=F)
#df1$Type=factor(df1$Type,levels=c("CTCF_U","U_specific","CTCF_I","I_specific"))
#png('WT.KO.LG.distance.png',width=6*500,height=7*500,res=500)
pdf(paste(prefix,'WT.KO.LG.distance.pdf',sep='.'))
#p=ggplot(df1,aes(x=dist,group=Type))+geom_density(aes(color=Type),size=1.5)+
#p=ggplot(df1,aes(x=Type,y=dist))+
p=ggplot(data,aes(x=Type,y=dist))+
geom_violin(aes(fill=Type),position = "dodge",size=.5)+
geom_boxplot(fill="white",width=.1,outlier.shape=NA,size=.5)+
scale_x_discrete(breaks=c("Lost","Retain","Gain"),labels=xlabel)+
scale_fill_manual(values=c("blue","#D0D0D0","red"))+
xlab("")+
ylab('Loop Size (kb)')+
#ggtitle('Loop size distribution')+
theme_bw()+
theme(legend.title = element_blank(),
        legend.text = element_text(size=14,face="bold"),
        legend.position = 'right') +
theme(axis.title=element_text(size=18,face="bold"),
        axis.text = element_text(color='black', size=10,face="bold"),
        axis.text.x = element_text(angle = 0,size=10,face='bold'),
        axis.text.y=element_text(angle=0,size=13,face='bold'),
        plot.title = element_text(size=20, face='bold')) +
theme(axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.border=element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(hjust = 0.5,face="bold"))
print(p);dev.off()

