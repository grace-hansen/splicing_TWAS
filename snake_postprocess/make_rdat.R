#!/usr/bin/Rscript
library(data.table)
library(tidyverse)
library(optparse)
library(intervals)
library(leafcutter)
library(foreach)
library(edgeR)

arguments <- parse_args(OptionParser(usage = "make_allreads_file.R <trait> <tissue> <offset> <data_source>",option_list=list()),
                        positional_arguments = 4)
opt=arguments$opt
trait=arguments$args[1]
tissue=arguments$args[2]
offset=as.numeric(arguments$args[3])
data_source=arguments$args[4]

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
  if (data_source=="GTEx_v7") {
    vals<-vals[4:nrow(vals),]
  } else if (data_source == "CMC" || data_so=="GTEx_v8")  {
    vals<-vals[2:nrow(vals),]
  }
  subject_IDs<-fread(paste(tissue,"_subject_IDs",sep=""),header=FALSE)
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
  allele_counts<-as.data.frame(t(allele_counts),stringsAsFactors = FALSE) %>% mutate(V2=as.numeric(V2)) %>% arrange(.,V2)
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
  genos<-genos[complete.cases(genos),]
  
  write.table(genos_untyped,paste("results/posthoc/test_diff_",intron,"_",rsid,".txt",sep=''),sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
  system(paste("chmod a+rw results/posthoc/test_diff_",intron,"_",rsid,".txt",sep=''))
  return(list(genos,genos_untyped))
}

################ Do the things #####################
#Get the first gene, use it to make rdat file
setwd(paste("/project2/nobrega/grace/splicing/",trait,"/",tissue,sep=''))
cmd<-paste("head -2 /project2/nobrega/grace/splicing/",trait,"/",tissue,"/results/",data_source,".all.dat.top.gene | tail -1",sep='')
out<-system(cmd,intern=TRUE)
intron<-strsplit(out,'\t')[[1]][3]
rsid<-strsplit(out,'\t')[[1]][8]

prep_plink_genos(intron)
genos<-make_leafcutter_genos(intron,rsid)
if (data_source=="GTEx_v7" || data_source=="CMC") {
  gencode<-"/project2/nobrega/grace/gencode.v19.annotation_hg19.txt.gz"
} else {
  gencode<-"/project2/nobrega/grace/gencode.v29.annotation_hg38.txt.gz"
}
cmd=paste("~/bin/leafcutter/scripts/leafcutter_ds.R -i 0 -g 0 -c 0 --num_threads 4 --min_samples_per_intron 1 --exon_file=",gencode," -o results/posthoc/",
          intron,"_",rsid," ",tissue,".leafcutter_perind.numers.gz results/posthoc/test_diff_",intron,"_",rsid,".txt",sep='')
system(cmd)

cmd=paste("~/bin/leafcutter/leafviz/prepare_results.R -m results/posthoc/test_diff_",
          intron,"_",rsid,".txt -f 1 -o results/posthoc/reads_dat/",tissue,".Rdat ",tissue,
          ".leafcutter_perind.numers.gz results/posthoc/",intron,"_",rsid,
          "_cluster_significance.txt results/posthoc/",intron,"_",rsid,
          "_effect_sizes.txt ~/bin/leafcutter/hg19/hg19",sep='')
print(cmd)
system(cmd)

cmd=(paste("rm results/posthoc/test_diff_",intron,"_",rsid,".txt results/posthoc/",intron,"_",rsid,
           "_allreads_cluster_significance.txt results/posthoc/",intron,"_",rsid,
           "_allreads_effect_sizes.txt",sep=''))
system(cmd)

write("","results/posthoc/reads_dat/made_allreads_file.txt")