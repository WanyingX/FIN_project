#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
k=as.numeric(args[1])
maxs=as.numeric(args[2])
file=args[3]
library(ggplot2)
library(MASS)
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

df=read.table(file,sep='\t',stringsAsFactors=F,header=T)
colnames(df)[1:2]=c('WT','KO')

df$density=get_density(df$WT,df$KO,n=2^8)
col=c(colorRampPalette(c("#1235e0","#ABDDA4"))(3),colorRampPalette(c("#E6F598","red"))(47))
fill_range=seq(min(df$density),max(df$density), len=50)
outfile=paste(file,'.noframe.png',sep='.')
#maxs=16
png(outfile,width=7*500,height=7*500,res=500)
#png(paste('scatter.liver.top',i,'png',sep='.'),width=7g*500,height=7*500,res=500)
p=ggplot()+
geom_point(df,mapping=aes(WT,KO,color=density),size=1)+
scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
coord_cartesian(xlim = c(0,maxs),ylim=c(0,maxs*k))+
scale_color_gradientn(colours = col,breaks=fill_range,space = "Lab")+
theme_bw()+
theme(  legend.title = element_blank(),
        legend.text = element_text(size=12,face="bold"),
        legend.position = 'none') +
theme(  axis.title=element_blank(),
        axis.text = element_blank(),
        plot.title = element_blank(),
        axis.ticks = element_blank())+
theme(axis.line = element_blank(),
#        panel.border=element_rect(colour = "black", fill=NA, size=1),
        panel.border=element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_blank())+
geom_abline(slope = k,intercept=0,color="grey60",size=1.5)
print(p);dev.off()
