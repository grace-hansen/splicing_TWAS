#!/usr/bin/Rscript
library(data.table)
library(stringr)
library(tidyverse)
library(optparse)
library(ggrepel)
library(intervals)
library(leafcutter)
library(foreach)
library(edgeR)

arguments <- parse_args(OptionParser(usage = "quantify_splicing_allreads.R <gene> <trait> <tissue>",option_list=list()),
                        positional_arguments = 4)
opt=arguments$opt
trait=arguments$args[1]
tissue=arguments$args[2]
gene=arguments$args[3]
data_source=arguments$args[4]
pony_colors<-fread("/project2/nobrega/grace/pony_palette")

offset<-50000

#For a given gene:
  #Identify most significant intron in TWAS, and the cluster associated with that intron
  #Identify the rsid with the highest weight in the model of the regulation of that intron
  #Run leafcutter_ds using that rsid as a genotype
  #Run prepare_results, with no thresholds for minimum number of samples with reads or minimum number of reads for intron
  #Using prepare results's output Rdat and genotype (in 3 levels, e.g. AA, AG, GG), quantify:
    #Number of reads in annotated vs. unannotated introns in the cluster of interest by genotype
    #Number of reads in annotated vs. unannotated introns in the gene by genotype
  #Plot all reads in the original cluster of interest by genotype (it'll be messy)

#Make plink from cluster
prep_plink_genos<-function(intron) {
  clustID=strsplit(intron,'-')[[1]][1]
  chrom=strsplit(intron,'-')[[1]][2]
  start=as.numeric(strsplit(intron,'-')[[1]][3])-offset
  stop=as.numeric(strsplit(intron,'-')[[1]][4])+offset
  
  cmd=paste("zcat /project2/nobrega/grace/splicing/",trait,"/",tissue,"/",tissue,".leafcutter_perind.decs.",
            chrom,".gz | grep ",chrom,":",as.numeric(strsplit(intron,'-')[[1]][3]),
            ":",as.numeric(strsplit(intron,'-')[[1]][4]),":",clustID,sep='')
  vals<-as.data.frame(strsplit(system(cmd,intern=TRUE),'[\t]'),stringsAsFactors = FALSE)
  if (data_source=="GTEx") {
    vals<-vals[4:nrow(vals),]
  } else {
    vals<-vals[3:nrow(vals),]
  }
  subject_IDs<-fread(paste("/project2/nobrega/grace/splicing/",trait,"/",tissue,"/",tissue,"_subject_IDs",sep=""),header=FALSE)
  pheno<-cbind(subject_IDs,vals)
  write.table(pheno,paste(intron,".pheno",sep=''),sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
  phenofile=paste(intron,".pheno",sep='')
  
  cmd=paste("plink --bfile /project2/nobrega/grace/genos/",data_source,"_plink/chr",chrom,"_rsids --pheno ",phenofile,
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
  #Select only noisy IDs
  #noisy_IDs<-scan(paste("/project2/nobrega/grace/noisy/",tissue,"_noisyIDs",sep=''),what="character",sep='\n')
  #genos_untyped<-genos_untyped[genos_untyped$V1 %in% noisy_IDs]
  genos<-as.data.frame(genos)
  #genos<-genos[genos$V1 %in% noisy_IDs,]
  return(list(genos,genos_untyped))
}

count_genos<-function(genos_typed) {
#Count number of reads by genotype and verdict
  reads_by_geno<-as.data.frame(matrix(nrow=3,ncol=3))
  geno_counts<-table(genos_typed$genos)[order(table(genos_typed$genos))]
  geno_counts<-geno_counts[(length(geno_counts)-2):length(geno_counts)]
  if(length(names(geno_counts))==2) { #E.g. if there are no minor heterozygotes
    maj_allele<-names(table(c(strsplit(names(geno_counts),'/')[[1]],strsplit(names(geno_counts),'/')[[2]]))[1])
    min_allele<-names(table(c(strsplit(names(geno_counts),'/')[[1]],strsplit(names(geno_counts),'/')[[2]]))[2])
    if (min_allele %in% strsplit(names(geno_counts),'/')[[1]]==TRUE) { #Get hetero naming convention as it appears in the data, e.g. "G/A" or "A/G"
      hetero<-names(geno_counts)[1]
    } else {
      hetero<-names(geno_counts[2])
    }
    major_homo<-paste(maj_allele,maj_allele,sep='/')
    minor_homo<-paste(min_allele,min_allele,sep='/')
  } else {
    if (strsplit(names(geno_counts[3]),"/")[[1]][1]!=strsplit(names(geno_counts[3]),"/")[[1]][2]) { #Heterozygous is most common genotype 
      minor_homo<-names(geno_counts[1])
      major_homo<-names(geno_counts[2])
      hetero<-names(geno_counts[3])
    } else {
      minor_homo<-names(geno_counts[1])
      major_homo<-names(geno_counts[3])
      hetero<-names(geno_counts[2])
    }
  }
  genos_typed$dosage[genos_typed$genos==major_homo]<-0
  genos_typed$dosage[genos_typed$genos==hetero]<-1
  genos_typed$dosage[genos_typed$genos==minor_homo]<-2
  return(genos_typed)
}

noisy_splicing_boxplot<-function(gene,dat,intron,rsid) {
  dat<-dat[dat$genos!='0/0',]
  G<-ggplot(dat,aes(x=genos,y=prop,fill=genos))+
    geom_boxplot(width=0.5,outlier.shape=NA)+
    geom_jitter(shape=16,size=1.5,width=0.1,height=0)+
    ylab("Proportion unannotated reads")+
    xlab(paste(rsid," genotype",sep=''))+
    scale_fill_manual(values=c(rgb(pony_colors[4,1:3]),rgb(pony_colors[3,1:3]),rgb(pony_colors[5,1:3])))+
    guides(fill=FALSE)+
    theme_minimal()+
    ggtitle(gene)+
    theme(axis.title.x=element_text(size=15),
          axis.text.x=element_text(size=15),
          axis.title.y=element_text(size=15),
          axis.text.y=element_text(size=15))
  ggsave(paste("boxplots/",gene,"_noisy_boxplot.pdf",sep=''),width=4,height=4)
}

quantify_splicing<-function(gene,intron,rsid,genos_typed,boxplot=TRUE) {
  #load(paste("/project2/nobrega/grace/noisy/reads_dat/",tissue,"_allreads.Rdat",sep=""))
  load(paste("/project2/nobrega/grace/splicing/",trait,"/",tissue,"/results/posthoc/reads_dat/",tissue,".Rdat",sep=""))
  out<-data.frame(matrix(ncol=4,nrow=1))
  ########Cluster-level comparison of number of annotated vs. unannotated reads##########

  chrom<-strsplit(intron,'-')[[1]][2]
  intron_start<-strsplit(intron,"-")[[1]][3]
  intron_end<-strsplit(intron,"-")[[1]][4]
  print(intron_start)
  print(intron_end)
  print(introns[1,])
  cluster<-introns$clusterID[introns$start==intron_start & introns$end==intron_end][1] ######New cluster numbering system in allreads file
  
  clust_counts<-counts[which(sapply(strsplit(rownames(counts),":"),'[[',4)==cluster),] #Get rows from counts from cluster of interest
  clust_counts<-t(clust_counts)
  clust_verdicts=matrix(nrow=ncol(clust_counts),ncol=2)
  
  #Find verdicts for all introns in cluster
  for (i in 1:ncol(clust_counts)) { 
    int<-colnames(clust_counts)[i]
    start<-strsplit(int,":")[[1]][2]
    stop<-strsplit(int,":")[[1]][3]
    clust<-gsub('.{2}$','',strsplit(int,":")[[1]][4])
    cmd=paste("zcat /home/gthansen/bin/leafcutter/hg19/hg19_all_introns.bed.gz | grep ",chrom," | grep ",start," | grep ",stop,sep="")
    v<-system(cmd,intern=TRUE)
    if (length(v)>0) {
      verdict="annotated"
    } else { verdict="unannotated" 
    }
    clust_verdicts[i,1:2]<-c(colnames(clust_counts)[i],verdict)
  }
  annotated_cols<-which(clust_verdicts[,2]=="annotated")
  unannotated_cols<-which(clust_verdicts[,2]!="annotated")
  
  clust_dat<-count_genos(genos_typed)
  clust_dat <- clust_dat[match(row.names(clust_counts), clust_dat$V1),]
  for (i in 1:nrow(clust_dat)){
    clust_dat$annotated[i]<-sum(clust_counts[i,annotated_cols]) #assumes clust_counts and clust_dat have same row order
    clust_dat$unannotated[i]<-sum(clust_counts[i,unannotated_cols])
    clust_dat$total[i]<-clust_dat$annotated[i]+clust_dat$unannotated[i]
    clust_dat$prop[i]<-clust_dat$unannotated[i]/sum(clust_dat$annotated[i]+clust_dat$unannotated[i])
  }

  #Is there a relationship between genotype and number of reads in cluster?
  tot_clust.lm<-lm(total~dosage,data=clust_dat)
  f<-summary(tot_clust.lm)$fstatistic
  p<-pf(f[1],f[2],f[3],lower.tail=F)
  out[1,1:2]<-c(f[1],p)
  colnames(out)[1:2]<-c("F_cluster_numreads","p_cluster_numreads")
  if (length(annotated_cols)!=0 & length(unannotated_cols)!=0) {
    if (boxplot==TRUE){
      plot_dat<-clust_dat[clust_dat$prop!="NaN",]
      noisy_splicing_boxplot(gene,plot_dat,intron,rsid)
    }
    #Is there a relationship between genotype and the proportion of reads that are unannotated?
    prop_unanno_clust.lm<-lm(prop~dosage,data=clust_dat)
    f<-summary(prop_unanno_clust.lm)$fstatistic
    p<-pf(f[1],f[2],f[3],lower.tail=F)
    out[1,3:4]<-c(f[1],p)
  }
  colnames(out)[3:4]<-c("F_prop_unannotated","p_prop_unannotated")
  write.table(out,paste("noisy/lm/",gene,"_",rsid,"_",cluster,"_lm.txt",sep=''),row.names=FALSE,quote=FALSE,sep='\t')
}

################ Do the things ###################

setwd(paste("/project2/nobrega/grace/splicing/",trait,"/",tissue,"/results/posthoc",sep=''))

#get weightfiles of most significant intron in each cluster
cmd=paste("awk -F '\t' '{ if ($1 == \"",gene,"\") print $2}' /project2/nobrega/grace/splicing/",trait,"/",tissue,"/results/v7.",data_source,".all.dat.top.gene | cut -f2 | cut -f5 -d'/'",sep="")
all_wgts<-system(cmd,intern=TRUE)
clusts<-sapply(strsplit(all_wgts,"-"),'[[',1) %>% unique()
weightfiles<-character()
for ( i in 1:length(clusts) ) {
  all_wgts[grep(clusts[i],all_wgts)[1]]
  weightfiles<-c(weightfiles,all_wgts[grep(clusts[i],all_wgts)[1]])
}

for ( i in 1:length(weightfiles)) {
  #Make gene plots and diff splicing plots for the first intron, make diff splicing plots on all other introns
  intron<-strsplit(weightfiles[i],'[.]')[[1]][1]
  
  #Get rsid: coloc'ed rsid if in coloc, top finemapped eQTL rsid if not in coloc
#  coloc_genes<-fread("coloc/colocalizing_genes_rsids",header=FALSE)
#  if (intron %in% coloc_genes$V2) {
#    rsid_dat<-fread(paste("coloc/",intron,"_coloc_full",sep=''))
#    rsid<-rsid_dat$snp[which.max(rsid_dat$SNP.PP.H4)]
#  } else if (length(which(list.files("coloc/")==paste("sQTL_",intron,"_coloc_finemap",sep='')))>0) {
#    rsid_dat<-fread(paste("coloc/sQTL_",intron,"_coloc_finemap",sep=''))
#    rsid<-rsid_dat$snp[which.max(rsid_dat$SNP.PP)] 
#  } else {
#    cmd=paste("awk -F '\t' '$2==\"",intron,"\" {print $9}' ../v7.",data_source,".all.dat.top",sep='')
#    rsid=system(cmd,intern=TRUE)
#  }
#  if (rsid=='null') {
#    cmd=paste("awk -F '\t' '$2==\"./",tissue,"/WEIGHTS//",intron,".wgt.RDat\" {print $10}' ../v7.",data_source,".all.dat.gene",sep='')
#    rsid=system(cmd,intern=TRUE)
#  }
cmd=paste("awk -F '\t' '$2==\"./",tissue,"/WEIGHTS//",intron,".wgt.RDat\" {print $10}' ../v7.",data_source,".all.dat.gene",sep='')
rsid=system(cmd,intern=TRUE)
   
  prep_plink_genos(intron)
  genos<-make_leafcutter_genos(intron,rsid)
  genos_typed<-genos[[1]]
  #Perform linear modelling  
  tryCatch({
    quantify_splicing(gene,intron,rsid,genos_typed)
  },error=function(e){})

  #Remove extra files
  system(paste("rm *",intron,".*",sep=''))
  write(paste(gene," complete"),paste("noisy/lm/",gene,"_status.txt",sep=''))
  write(paste(gene," complete"),paste("boxplots/",gene,"_status.txt",sep=''))
}
