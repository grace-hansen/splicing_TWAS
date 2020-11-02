#!/usr/bin/Rscript
args <- commandArgs(trailingOnly = TRUE)
library(data.table)
library(tidyverse)
library(gap)
library(gridExtra)
library(yarrr)
trait<-args[1]
tissue<-args[2]
setwd(paste("~/midway/splicing/",trait,"/",tissue,"/results/posthoc/noisy/",sep=''))
genes_clusts<-fread("../../v7.GTEx.all.dat.top.gene",sep='\t',drop=seq(3,20)) #Get list of all sig genes and all sig clusts for each gene
genes_clusts$clust<-sapply(strsplit(sapply(strsplit(genes_clusts$FILE,'/'),'[[',5),'-'),'[[',1)
genes_clusts$FILE<-NULL
genes_clusts<-genes_clusts %>% distinct() %>% na.omit()

#Collect p values and F statistics for all genes in a single data frame
lm_dat<-data.frame(matrix(nrow=0,ncol=6))

for (i in 1:nrow(genes_clusts)) {
  #print(paste(i,genes_clusts$GENE[i],genes_clusts$clust))
  print(i)
  #Read in lm, try to add to lm_dat
  lm_files<-list.files("lm/",pattern="*lm.txt")[which(grepl(paste(genes_clusts$GENE[i],"_",sep=''),list.files("lm/",pattern="*lm.txt")))]
  if (length(lm_files) != 0) {
    for (j in 1:length(lm_files)) {
      lm<-fread(paste("lm/",lm_files[j],sep=''))
      lm_dat<-rbind(lm_dat,cbind(genes_clusts[i,1],lm_files[i],lm),use.names=FALSE) 
    }
  }
}
colnames(lm_dat)<-c("gene","clust","F_cluster_numreads","p_cluster_numreads",
                    "F_cluster_unannotated","p_cluster_unannotated")

write.table(lm_dat,"lm_no_covars.txt",row.names=FALSE,quote=FALSE,sep='\t')

qqplots<-function(dat,name) {
  pony<-piratepal(palette="pony")
  
  int<-lm(quantile(-log10(dat$p_cluster_unannotated),c(0.25,0.5),na.rm=TRUE)~c(0.25,0.5))$coeff[[1]]
  coef<-lm(quantile(-log10(dat$p_cluster_unannotated),c(0.25,0.5),na.rm=TRUE)~c(0.25,0.5))$coeff[[2]]
  A<-ggplot(dat,aes(sample=-log10(p_cluster_unannotated)))+
    stat_qq(color=pony[1],distribution=stats::qunif)+
    geom_abline(intercept=int,slope=coef)+
    ggtitle("Proportion of unannotated reads in cluster")+
    scale_y_continuous(name="-log10(p)",limits=c(0,40))+
    scale_x_continuous(name="Quantiles, uniform distribution")+
    theme_minimal()+
    theme(legend.position="none")
  
  int<-lm(quantile(-log10(dat$p_cluster_numreads),c(0.25,0.5),na.rm=TRUE)~c(0.25,0.5))$coeff[[1]] #Having to calculate qqline by hand, b/c stat_qq_line is buggy with scale limits
  coef<-lm(quantile(-log10(dat$p_cluster_numreads),c(0.25,0.5),na.rm=TRUE)~c(0.25,0.5))$coeff[[2]]
  B<-ggplot(dat,aes(sample=-log10(p_cluster_numreads)))+
    stat_qq(aes(color=pony[8]),distribution=stats::qunif)+
    geom_abline(intercept=int,slope=coef)+
    ggtitle("Number of reads in cluster")+
    scale_y_continuous(name="-log10(p)",limits=c(0,40))+
    scale_x_continuous(name="Quantiles, uniform distribution")+
    theme_minimal()+
    theme(legend.position="none")

  g<-arrangeGrob(A,B,nrow=1)
  
  ggsave(paste("qqplots_",name,".pdf",sep=''),width=8,height=4,g)
}

##Make qqplot, no PCAs
qqplots(lm_dat,paste("splicing_",trait,"_",tissue,sep=''))