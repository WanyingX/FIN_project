#!/bin/bash
# check /mnt/rds/genetics01/JinLab/xww/0326_data/09-26-2023-pixel-gene/10-10-2023

#bedtools intersect -wa -wb -a ~/HiC/references/mm10.HindIII.anchors_5kb.bed -b ~/RNA-seq/mm10.refFlat.TSS.bed | awk '{print $9,$4}' OFS='\t' | sed 's/\t/:/g' | sort | uniq | sed 's/:/\t/g' > HindIII.gene.2.anchor

#bedtools intersect -wa -wb -a ~/HiCorr/references/mm10_DpnII_anchors_avg.bed -b ~/RNA-seq/mm10.refFlat.TSS.bed | awk '{print $9,$4}' OFS='\t' | sed 's/\t/:/g' | sort | uniq | sed 's/:/\t/g' > DpnII.gene.2.anchor

#bedtools intersect -wa -wb -a ~/micro_C/ref/mm10.genome.5kb.bed -b ~/RNA-seq/mm10.refFlat.TSS.bed | awk '{print $9,$4}' OFS='\t' | sed 's/\t/:/g' | sort | uniq | sed 's/:/\t/g' > microc.gene.2.anchor

echo FIN `cat FIN.ctcf.peak.bed | awk '{print $1,$2-100000,$3+100000}' OFS='\t' | awk '{if($2>0) print $0;else print $1,0,$3}' OFS='\t' | cut -f1-3 | sort -k1,1 -k2n | bedtools merge  | awk '{s+=$3-$2}END{print s}'` >> cover.stat

echo TAD `cat mm10.TAD.boundary | awk '{print $1,$2-100000,$3+100000}' OFS='\t' | awk '{if($2>0) print $0;else print $1,0,$3}' OFS='\t' | cut -f1-3 | sort -k1,1 -k2n | bedtools merge  | awk '{s+=$3-$2}END{print s}'` >> cover.stat

echo pixel `cat recurrent.gainloop.anchor  | cut -f1-3 | sort -k1,1 -k2n | bedtools merge  | awk '{s+=$3-$2}END{print s}'` >> cover.stat




direct=$1


for file in `ls *.DEGS.table`;do

echo $file FIN gene `cat FIN.ctcf.peak.bed | awk '{print $1,$2-100000,$3+100000}' OFS='\t' | awk '{if($2>0) print $0;else print $1,0,$3}' OFS='\t' | bedtools intersect -wa -wb -a - -b $file | awk '{print $(NF-1),$NF}' OFS=':' | sort | uniq | sed 's/:/\t/g' | grep -w $direct |wc -l` >> FIN.ctcf.$direct.stat
#echo $file FIN cover `cat FIN.ctcf.peak.bed | awk '{print $1,$2-100000,$3+100000}' OFS='\t' | awk '{if($2>0) print $0;else print $1,0,$3}' OFS='\t' | bedtools intersect -wa -wb -a - -b $file | cut -f1-3 | sort -k1,1 -k2n | bedtools merge | awk '{s+=$3-$2}END{print s}'` >> FIN.ctcf.stat
echo $file FIN cover 161289079 >> FIN.ctcf.$direct.stat
done &


for file in `ls *.DEGS.table`;do

echo $file TAD gene `cat mm10.TAD.boundary | awk '{print $1,$2-100000,$3+100000}' OFS='\t' | awk '{if($2>0) print $0;else print $1,0,$3}' OFS='\t' | bedtools intersect -wa -wb -a - -b $file | awk '{print $(NF-1),$NF}' OFS=':' | sort | uniq | sed 's/:/\t/g' | grep -w $direct |wc -l` >> TAD.ctcf.$direct.stat
#echo $file FIN cover `cat mm10.TAD.boundary | awk '{print $1,$2-100000,$3+100000}' OFS='\t' | awk '{if($2>0) print $0;else print $1,0,$3}' OFS='\t' | bedtools intersect -wa -wb -a - -b $file | cut -f1-3 | sort -k1,1 -k2n | bedtools merge | awk '{s+=$3-$2}END{print s}'` >> TAD.ctcf.stat
echo $file TAD cover 713600000 >> TAD.ctcf.$direct.stat
done &


for file in `ls *.DEGS.table`;do

echo $file gain_pixel gene `cat recurrent.gainloop.anchor | bedtools intersect -wa -wb -a - -b $file | awk '{print $(NF-1),$NF}' OFS=':' | sort | uniq | sed 's/:/\t/g' | grep -w $direct | wc -l` >> gainpixel.$direct.stat
#echo $file FIN cover `cat recurrent.gainloop.anchor| bedtools intersect -wa -wb -a - -b $file | cut -f1-3 | sort -k1,1 -k2n | bedtools merge | awk '{s+=$3-$2}END{print s}'` >> pixel.stat

echo $file gain_pixel cover 16016796 >> gainpixel.$direct.stat

done &


for file in `ls *.DEGS.table`;do

echo $file lost_pixel gene `cat nora.lostloop.anchor | bedtools intersect -wa -wb -a - -b $file | awk '{print $(NF-1),$NF}' OFS=':' | sort | uniq | sed 's/:/\t/g' | grep -w $direct | wc -l` >> lostpixel.$direct.stat
#echo $file FIN cover `cat recurrent.gainloop.anchor| bedtools intersect -wa -wb -a - -b $file | cut -f1-3 | sort -k1,1 -k2n | bedtools merge | awk '{s+=$3-$2}END{print s}'` >> pixel.stat

echo $file lost_pixel cover 277531307 >> lostpixel.$direct.stat

done &

for file in `ls *.DEGS.table`;do

echo $file common_pixel gene `cat nora.commonloop.anchor | bedtools intersect -wa -wb -a - -b $file | awk '{print $(NF-1),$NF}' OFS=':' | sort | uniq | sed 's/:/\t/g' | grep -w $direct | wc -l` >> commonpixel.$direct.stat
#echo $file FIN cover `cat recurrent.gainloop.anchor| bedtools intersect -wa -wb -a - -b $file | cut -f1-3 | sort -k1,1 -k2n | bedtools merge | awk '{s+=$3-$2}END{print s}'` >> pixel.stat

echo $file common_pixel cover 99192287 >> commonpixel.$direct.stat

done &


for d in down up;do

cat *$d*.stat  | grep nora | sed 's/\./ /g' > nora.$d.stat &
cat *$d*.stat  | grep kubo | sed 's/\./ /g' > kubo.$d.stat &
cat *$d*.stat | grep hsieh.mNET | sed 's/\./ /g' > hsieh.mNET.$d.stat &
cat *$d*.stat | grep hsieh.RNA | sed 's/\./ /g' > hsieh.RNA.$d.stat &
done

#filelist=c()
direct=c("up","down")

name=c("nora","kubo","hsieh.mNET","hsieh.RNA")
filelist=paste(name,"stat",sep='.')

for(dd in 1:length(direct)){

filelist=paste(name,direct[dd],"stat",sep='.')

pdf(paste(direct[dd],".gene.coverage.pdf",sep=''),width=6,height=8)

nf=layout(matrix(1:12,4,3,byrow=TRUE),width=rep(3,12),height=rep(5,12),TRUE)


par(mar=c(2,2,2,2))

for(i in 1:length(name)){

	df=read.table(filelist[i],sep=' ',stringsAsFactors=F)
	if(ncol(df)==9){
	df[,1]=NULL
	}
	gene=df[df[,7]=="gene",]
	cover=df[df[,7]=="cover",]
	colnames(gene)[8]="gene"
	colnames(cover)[8]="coverage"
	cover$coverage=cover$coverage/1000000
	time=names(table(gene[,3]))
	out=data.frame(time=gene[,3],group=gene[,6],gene=gene$gene,coverage=cover$coverage*1000000,study=name[i],val=gene$gene/cover$coverage)

	write.table(out,paste(name[i],direct[dd],"table",sep='.'),sep='\t',quote=F,row.names=F)

	data=data.frame(time=gene[,3],group=gene[,6],gene=gene$gene,coverage=cover$coverage)
	data$val=data$gene/data$coverage
	data$group=factor(data$group,levels=c("TAD","FIN","lost_pixel","common_pixel","gain_pixel"))
	for(tt in time){
		m = data[data$time==tt,]
		rownames(m)=m$group
#		barplot(m$val~m$group,col=c("orange","darkgreen","red"),main=paste(name[i],tt,sep='-'))
		maxs=max(m$val)+0.2
#                barplot(m$val~m$group,col=c("#00AEBA","#E5B721","#EF5123"),border=c("#00AEBA","#E5B721","#EF5123"),main=paste(name[i],tt,sep='-'),width=rep(0.2,3),ylim=c(0,maxs))
		barplot(m$val~m$group,col=c("#00AEBA","#E5B721","#533cc7","#959c97","#EF5123"),border=c("#00AEBA","#E5B721","#533cc7","#959c97","#EF5123"),main=paste(name[i],tt,sep='-'),width=rep(0.2,3),ylim=c(0,maxs))
	}

}

dev.off()
}



mat=data.frame()

direct=c("up","down")
name=c("nora","kubo","hsieh.mNET","hsieh.RNA")
#filelist=paste(name,"stat",sep='.')
for(dd in 1:length(direct)){
	filelist=paste(name,direct[dd],"table",sep='.')
	for(i in 1:length(name)){

        df=read.table(filelist[i],sep='\t',stringsAsFactors=F,header=T)
	df$direct=direct[dd]
  	mat=rbind(mat,df)
 }
}

b1=c("1day","3hrs")
b2=c("2days","12hrs")
b3=c("4days","24hrs")
data=mat[mat$time%in%b1,]
data$val=(data$gene/data$coverage)*1000000
#        data=data[data$group!="pixel",]
data$group=factor(data$group,levels=c("TAD","FIN","lost_pixel","common_pixel","gain_pixel"))
data$study=factor(data$study,levels=name)
library(ggplot2)

for(dd in direct){

pdf(paste(dd,"early.pdf",sep='.'),7,7)

m=data[data$direct==dd,]

p=ggplot(m,mapping=aes(x=study,y=val,fill=group))+
#geom_point(df,mapping=aes(dist,residual),color='grey',size=.5)+
geom_bar(stat="identity", width=.6,position="dodge")+
#facet_wrap(~type)+
#xlab('Types of loops')+
ylab('Enrichment score')+
#ggtitle('Enrichment score of histone marker')+
#geom_text(data=anno,aes(x,y,label=label),size=5,color="black")+
xlab("")+
#scale_fill_manual(values = c("#00AEBA","#E5B721","#EF5123"))+
scale_fill_manual(values = c("#00AEBA","#E5B721","#533cc7","#959c97","#EF5123"))+

theme_bw()+
theme(legend.title = element_blank(),
        legend.text = element_text(size=12,face="bold"),
        legend.position = 'bottom') +
theme(  #axis.title=element_blank(),
        axis.title=element_text(color='black', size=10,face="bold"),
        axis.text = element_text(color='black', size=7,face="bold"),
        #axis.text = element_blank(),
        axis.text.x = element_text(angle = 0,size=13,face='bold'),
        axis.text.y=element_text(angle=0,size=13,face='bold'),
        plot.title = element_text(size=10, face='bold')) +
theme(axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.border=element_rect(colour = "black", fill=NA, size=1),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(hjust = 0.5,face="bold"))

print(p);dev.off()
}
