#!/bin/env Rscript
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)

name=args[1]

markerep=read.table(args[2],sep='\t',stringsAsFactors=F)

col=paste("#",c('9197b5',"d8b23c",'e21f26',"2b2f84","82b297","bdc3c5"),sep='')
names(col)=c('CTCF-loop',"CTCF_EP",'EP_only','K27me3','K36me3',"Other")

inputfile=c('lost.loops','common.loops','gain.loops')
names(inputfile)=c('lost','common','gain')
file=c("CTCF.relevant.loops","CTCF_CTCF",'merged.CTCF','merged.EP','merged.H3K27me3','merged.H3K36me3')
names(file)=c("All","CTCF-CTCF",'CTCF-marker','EP','K27me3','K36me3')

datalist=list()

all=read.table("CTCF.relevant.loops",sep='\t',stringsAsFactors=F)
ctcf=read.table('merged.CTCF',sep='\t',stringsAsFactors=F)
only=read.table("CTCF_CTCF",sep='\t',stringsAsFactors=F)
ep=read.table('merged.EP',sep='\t',stringsAsFactors=F)


rownames(all)=paste(all[,1],all[,2],sep=':')
rownames(only)=paste(only[,1],only[,2],sep=':')
rownames(ctcf)=paste(ctcf[,1],ctcf[,2],sep=':')
rownames(ep)=paste(ep[,1],ep[,2],sep=':')

datalist[["CTCF_EP"]]=intersect(rownames(all),rownames(ep))
#datalist[['CTCF_only']]=setdiff(rownames(ctcf),intersect(rownames(ctcf),rownames(ep)))
datalist[['EP_only']]=setdiff(rownames(ep),intersect(rownames(all),rownames(ep)))

all$group="CTCF-unknown"
all[rownames(ctcf),"group"]="CTCF-other"
all[rownames(only),"group"]="CTCF-CTCF"
all[setdiff(datalist[["CTCF_EP"]],intersect(rownames(only),datalist[["CTCF_EP"]])),"group"]="CTCF_EP"

datalist[["CTCF-CTCF"]]=rownames(all)[all$group=="CTCF-CTCF"]
datalist[["CTCF-other"]]=rownames(all)[all$group=="CTCF-other"]
datalist[["CTCF-unknown"]]=rownames(all)[all$group=="CTCF-unknown"]

k27=read.table("merged.H3K27me3",sep='\t',,stringsAsFactors=F)
k36=read.table("merged.H3K36me3",sep='\t',,stringsAsFactors=F)

k27=paste(k27[,1],k27[,2],sep=':')
k36=paste(k36[,1],k36[,2],sep=':')


k27=setdiff(k27,intersect(k27,union(rownames(all),datalist[['EP_only']])))
k36=setdiff(k36,intersect(k36,union(rownames(all),union(datalist[['EP_only']],k27))))

datalist[["K27me3"]]=k27
datalist[["K36me3"]]=k36

orders=c('CTCF-CTCF',"CTCF_EP",'EP_only',"CTCF-other","CTCF-unknown",'K27me3','K36me3')

out=data.frame()
for(i in names(inputfile)){


df=read.table(inputfile[i],sep='\t',stringsAsFactors=F,row.names=1)

df$Type="Other"
for(j in orders){

df1=datalist[[j]]
df[intersect(rownames(df)[which(df$Type=="Other")],df1),'Type']=j

}

temp=data.frame(matrix(unlist(strsplit(rownames(df),":")),nrow(df),2,byrow=T))
rownames(temp)=paste(temp[,1],temp[,2],sep=':')
temp$Type=df$Type
other=temp[which(temp$Type=="Other"),]

#other[(other[,1]%in%markerep[,4])||(other[,2]%in%markerep[,4]),"Type"]="EP_only"
other[other[,1]%in%markerep[,4],"Type"]="EP_only"
other[other[,2]%in%markerep[,4],"Type"]="EP_only"

df[rownames(other),"Type"]=other$Type
ctcfloop=c('CTCF-CTCF',"CTCF-other","CTCF-unknown")

df[df$Type%in%ctcfloop,"Type"]="CTCF-loop"

piedata=data.frame(Type=c('CTCF-loop',"CTCF_EP",'EP_only','K27me3','K36me3',"Other"),Freq=0)
rownames(piedata)=piedata$Type
tt=data.frame(table(df$Type))
rownames(tt)=tt[,1]
piedata[rownames(tt),'Freq']=tt$Freq
#piedata=data.frame(table(df$Type))

piedata$percent=piedata$Freq/nrow(df)
names(piedata)[1]='Type'
#print(piedata)
piedata$Type=factor(piedata$Type,levels=c('CTCF-loop',"CTCF_EP",'EP_only','K27me3','K36me3',"Other"))

piedata$group=i
out=rbind(out,piedata)
#pdf("lost.loops.pie.pdf")

outfile=paste(name,paste(i,'.pie.pdf',sep=''),sep='.')
pdf(outfile)
p=ggplot(piedata, aes("", percent, fill = Type)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0,direction=1)+
#  geom_text(aes(y = pos, label =label ),size=6,color = "white")+
geom_text(aes(label = paste(round(percent*100), "%",sep=''), x = 1.3),
            position = position_stack(vjust = 0.5))+
  scale_fill_manual(values = col) +
#scale_fill_manual(values = col[c(1,2,3,5,6)]) +
theme(  legend.title = element_blank(),
        legend.text = element_text(size=12,face="bold"),
        legend.position = 'none') +
  theme_void()
print(p);dev.off()


}

#write.table(out,'pie.stat',sep='\t',quote=F,row.names=F)

