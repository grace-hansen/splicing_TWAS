#!/usr/bin/Rscript
source("~/bin/leafcutter/leafcutter/R/make_gene_plot_bare.R")
library(data.table)
library(stringr)
library(tidyverse)
library(optparse)
library(ggrepel)
library(intervals)
library(leafcutter)
library(foreach)

arguments <- parse_args(OptionParser(usage = "plot_splicing.R <gene> <trait> <tissue> <data_source>",option_list=list()),
                        positional_arguments = 4)
opt=arguments$opt
gene=arguments$args[1]
trait=arguments$args[2]
tissue=arguments$args[3]
data_source=arguments$args[4]
offset<-50000

#Make plink from cluster
prep_plink_genos<-function(intron) {
  clustID=strsplit(intron,'-')[[1]][1]
  chrom=strsplit(intron,'-')[[1]][2]
  start=as.numeric(strsplit(intron,'-')[[1]][3])-offset
  stop=as.numeric(strsplit(intron,'-')[[1]][4])+offset
  
  cmd=paste("zcat /project2/nobrega/grace/splicing/",trait,"/",tissue,"/",tissue,"_perind.decs.",
            chrom,".gz | grep ",chrom,":",as.numeric(strsplit(intron,'-')[[1]][3]),
            ":",as.numeric(strsplit(intron,'-')[[1]][4]),":",clustID,sep='')
  vals<-as.data.frame(strsplit(system(cmd,intern=TRUE),'[\t]'),stringsAsFactors = FALSE)
  vals<-vals[2:nrow(vals),]
  subject_IDs<-fread(paste("/project2/nobrega/grace/splicing/",trait,"/",tissue,"/",tissue,"_subject_IDs",sep=""),header=FALSE)
  pheno<-cbind(subject_IDs,vals)
  write.table(pheno,paste(intron,".pheno",sep=''),sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
  phenofile=paste(intron,".pheno",sep='')
  
  cmd=paste("plink --bfile /project2/nobrega/grace/genos/",data_source,"/plink/chr",chrom,"_rsids --pheno ",phenofile,
            " --make-bed --out ",intron," --keep ",phenofile," --chr ",chrom,
            " --from-bp ",start," --to-bp ",stop,sep="")
  system(cmd)
  cmd=paste("plink --bfile ",intron," --recode tab --out ",intron,sep='')
  system(cmd)
}

#Read in SNP genotypes from plink
make_leafcutter_genos<-function(intron,rsid) {
  plink_rsids<-fread(paste(intron,".map",sep=''))
  line<-which(plink_rsids$V2==rsid)
  plink_genos<-fread(paste(intron,".ped",sep=''))
  IDs<-plink_genos[,1]
  plink_genos<-plink_genos[,7:ncol(plink_genos)]
  genos<-as.data.frame(plink_genos)[,line]
  genos<-gsub(" ", "/", genos)
  genos<-cbind(IDs,genos)
  allele_counts<-matrix(nrow=2,ncol=0)
  for (c in c("A","C","T","G")) {
    num<-sum(str_count(genos$genos,c))
    allele_counts<-cbind(allele_counts,as.data.frame(c(c,num),stringsAsFactors = FALSE))
  }
  allele_counts<-as.data.frame(t(allele_counts),stringsAsFactors = FALSE) %>% arrange(.,V2)
  major=allele_counts[4,1]
  minor=allele_counts[3,1]
  genos_untyped<-genos
  genos_untyped$genos<-NA
  
  if (paste(minor,minor,sep='/') %in% genos$genos) { #minor allele homozygotes are present
    for (i in 1:nrow(genos)) {
      if (genos$genos[i]==paste(major,'/',major,sep='') || genos$genos[i]==paste(major,'/',minor,sep='') ||  genos$genos[i]==paste(minor,'/',major,sep='')) {
        genos_untyped$genos[i]<-'major_carrier'
      } else if (genos$genos[i]==paste(minor,'/',minor,sep='')) {
        genos_untyped$genos[i]<-'minor_homo'
      }
    }
  } else {
    for (i in 1:nrow(genos)) {
      if (genos$genos[i]==paste(major,'/',major,sep='')) {
        genos_untyped$genos[i]<-'major_homo'
      } else if (genos$genos[i]==paste(minor,'/',major,sep='') || genos$genos[i]==paste(major,'/',minor,sep='')) {
        genos_untyped$genos[i]<-'hetero'
      }
    }
  }
  
  genos_untyped<-genos_untyped[complete.cases(genos_untyped),]
  write.table(genos_untyped,paste("test_diff_",intron,"_",rsid,".txt",sep=''),sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
  system(paste("chmod a+rw test_diff_",intron,"_",rsid,".txt",sep=''))
  return(genos_untyped)
}


#Plot results: two plots per gene(one with clusters, one without), one plot per significant cluster
plot_gene<-function(gene,intron) {  
  #Make gene plots
  G<-make_gene_plot_bare(gene,cluster_list=clusters,introns=introns,introns_to_plot=introns_to_plot,exons_table=exons_table)
  G_C<-make_gene_plot_bare(gene,cluster_list=clusters,show_clusters=TRUE,introns=introns,introns_to_plot=introns_to_plot,exons_table=exons_table)
  ggsave(paste("diff_splicing/",gene,"_isoform_plot.pdf",sep=""),plot=G,device="pdf",width=6,height=3)
  ggsave(paste("diff_splicing/",gene,"_isoform_plot_clusts.pdf",sep=""),plot=G_C,device="pdf",width=6,height=3)
}

#Plot results: two plots per gene(one with clusters, one without), one plot per significant cluster
plot_diff_splicing<-function(gene,weightfile,intron,rsid) {
  #Make cluster plots
  ###Get cluster ID in allreads file: won't be the same as clusterID from weightfile
  cluster<-strsplit(weightfile,'-')[[1]][1]
  start<-strsplit(weightfile,'-')[[1]][3]
  stop<-strsplit(strsplit(weightfile,'-')[[1]][4],'[.]')[[1]][1]
  C<-make_cluster_plot(cluster,exons_table=exons_table,meta=meta,cluster_ids=cluster_ids,counts=counts,introns=introns)
  ggsave(paste("diff_splicing/",gene,"_",intron,"_",rsid,"_diffsplicing_plot.pdf",sep=""),plot=C,device="pdf",width=6,height=6)
}

################ Do the things ###################

setwd(paste("/project2/nobrega/grace/splicing/",trait,"/",tissue,"/results/posthoc",sep=''))

#get introns from gene
cmd=paste("awk -F '\t' '{ if ($1 == \"",gene,"\") print $2}' /project2/nobrega/grace/splicing/",trait,"/",tissue,"/results/",data_source,".all.dat.top.gene | cut -f2 | cut -f5 -d'/'",sep="")
weightfiles<-system(cmd,intern=TRUE)

#Make gene plots and diff splicing plots for the first intron, make diff splicing plots on all other introns
for (i in 1:length(weightfiles)) {
  intron<-strsplit(weightfiles[i],'[.]')[[1]][1]

  #Get rsid: coloc'ed rsid if in coloc, top finemapped eQTL rsid if not in coloc
  coloc_genes<-fread("coloc/colocalizing_genes_rsids",header=FALSE)
  if (intron %in% coloc_genes$V2) {
    rsid_dat<-fread(paste("coloc/",intron,"_coloc_full",sep=''))
    rsid<-rsid_dat$snp[which.max(rsid_dat$SNP.PP.H4)]
  } else if (length(which(list.files("coloc/")==paste("sQTL_",intron,"_coloc_finemap",sep='')))>0) {
    rsid_dat<-fread(paste("coloc/sQTL_",intron,"_coloc_finemap",sep=''))
    rsid<-rsid_dat$snp[which.max(rsid_dat$SNP.PP)] 
  } else {
    cmd=paste("awk -F '\t' '$2==\"",intron,"\" {print $9}' ../",data_source,".all.dat.top",sep='')
    rsid=system(cmd,intern=TRUE)
  }
  if (rsid=='null') {
    cmd=paste("awk -F '\t' '$2==\"./",tissue,"/WEIGHTS/",intron,".wgt.RDat\" {print $10}' ../",data_source,".all.dat.gene",sep='')
    rsid=system(cmd,intern=TRUE)
  }
  
  cmd=paste("awk -F '\t' '$2==\"",intron,"\" {print $9}' ../",data_source,".all.dat.top",sep='')
  rsid=system(cmd,intern=TRUE)
  
  prep_plink_genos(intron)
  genos<-make_leafcutter_genos(intron,rsid)
  cmd<-paste("rm ",intron,"*",sep="")
  system(cmd)
  
  #Load Rdat
  load(paste("/project2/nobrega/grace/splicing/",trait,"/",tissue,"/results/posthoc/reads_dat/",tissue,".Rdat",sep=''))
  
  tryCatch({
    if (i==1) { 
      plot_gene(gene,intron)
    }
    plot_diff_splicing(gene,weightfiles[i],intron,rsid)
  },error=function(e){})
  
  cmd=paste("rm test_diff_",intron,"_",rsid,".txt",sep='')
  system(cmd)
}

