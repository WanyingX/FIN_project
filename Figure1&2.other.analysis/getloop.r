#!/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# This script is used to select specific loops. e.g. lost loops, gain loops and retain loops.
# 1. For lost and gain loops, we use 4 fold change approach. 
# 2. For retain loop, we use residuals, xintercept and yintecepet together to pick up retain loop.

k=as.numeric(args[1])
filename=args[2]
library(ggplot2)
library(MASS)
df=read.table(filename,sep='\t',stringsAsFactors=F,header=T,row.names=1)
a=df[order(-df$WT),]
#a=df[order(-df[,1]),]
df$residuals=df$KO-df$WT*k
#xintercept=min(a[1:100000,1])
xintercept=min(a[1:500000,1])
yintercept=xintercept*k
cutoff=k*xintercept
common=df[which(abs(df$residuals)<cutoff&df$WT>xintercept&df$KO>yintercept),]

df$fc=df$KO/df$WT
gain=df[df$fc>4,]
lost=df[df$fc<(1/4),]

write.table(lost,"lost.loops",sep='\t',quote=F,col.names=F)
write.table(gain,"gain.loops",sep='\t',quote=F,col.names=F)
write.table(common,'common.loops',sep='\t',quote=F)

