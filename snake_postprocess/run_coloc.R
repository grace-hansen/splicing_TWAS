#!/usr/bin/Rscript
args=commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("USAGE: run_coloc.R <trait> <tissue> <intron> <offset> <exp_N>\n", call.=FALSE)
}
library(coloc)
library(data.table)
library(tidyverse)
library(gridExtra)
source ("/project2/nobrega/grace/expression/scripts/coloc/claudia.R")
trait=args[1]
tissue=args[2]
intron=args[3]
window=args[4]
exp_N=args[5]

setwd(paste("/project2/nobrega/grace/splicing/",trait,"/",tissue,"/results/posthoc",sep=''))

#Get chrom
chrom<-strsplit(intron,'-')[[1]][2]

#Grab p values for both datasets
get_coloc_data<-function(dataset1_name,dataset2_name) { ##Names must be filenames in the same directory as wd
  tryCatch({
    dataset1<-fread(dataset1_name)
    dataset2<-fread(dataset2_name)
  
  ##Make sure both datasets have a p-value or beta
    get_P_or_B<-function(dataset) {
      if ("P" %in% colnames(dataset)) {
        test_type="P"
      } else {
        if ("B" %in% colnames(dataset)) {
          test_type="B"
        } else if ("Z" %in% colnames(dataset)) {
          dataset$P=2*pnorm(-abs(dataset$Z))
          test_type="P"
        } else {
          print("p value or beta cannot be found for one of the datasets, cannot run coloc")
        }
      }
      return(list(data=dataset,test_type=test_type))
    }

    #Find N for both datasets
    find_N<-function(dataset,dataset_name) {
      if ("N" %in% colnames(dataset)) {
        N<-dataset$N
      } else {
        N<-as.numeric(exp_N) #If there isn't an 'N' column in the data, there needs to be a [filename]_N file from which to read the N
      }
      return(N)
    }

    dset1<-get_P_or_B(dataset1)$data
    dset1_test_type<-get_P_or_B(dataset1)$test_type
    dset2<-get_P_or_B(dataset2)$data
    dset2_test_type<-get_P_or_B(dataset2)$test_type
    dset1_N<-find_N(dataset1,dataset1_name)
    dset2_N<-find_N(dataset2,dataset2_name)
  
    if (dset1_test_type==dset2_test_type) {
      return(list(dset1=dset1,dset2=dset2,dset1_N=dset1_N,dset2_N=dset2_N,test_type=dset1_test_type))
    } else {
      print("One dataset has p-values and one has beta values; this needs to be fixed before coloc can be run")
    }
  },
  error=function(e) {
    return(list(dset1=data.frame(),dset2=data.frame(),dset1_N='',dset2_N='',test_type=''))
  })
}

data<-get_coloc_data(paste("coloc/loci/",intron,"_GWAS_sumstats",sep=''),paste("coloc/loci/",intron,"_sQTL.txt",sep=''))
GWAS<-data$dset1
sQTL<-data$dset2
sQTL_N<-data$dset2_N
test_type<-data$test_type

if (nrow(GWAS)>0 && nrow(sQTL>0)) {
  #Window data to desired range
  start<-strsplit(intron,'-')[[1]][3]
  stop<-strsplit(intron,'-')[[1]][4]
  strand<-strsplit(data$dset2$intron_id[1],':')[[1]][5]
  if (strand=="+") {
    fiveprime=as.numeric(start)
  } else {
    fiveprime=as.numeric(stop)
  }
  
  sQTL$distance<-sQTL$distance
  GWAS$distance<-sQTL$distance
  sQTL$genom_distance<-sQTL$distance+fiveprime
  GWAS$genom_distance<-sQTL$distance+fiveprime
  
  GWAS<-GWAS[abs(GWAS$distance)<=window,]
  sQTL<-sQTL[abs(sQTL$distance)<=window,]
  
  ##Determine which snp is the likely causal variant in each dataset (we're assuming 1 here, which is probs wrong)
  MAF<-sQTL$MAF
  sQTL_finemap<-finemap.abf(dataset=list(pvalues=sQTL$P, N=sQTL_N, MAF=MAF, snp=sQTL$rsid,type="quant"))
  GWAS_finemap<-finemap.abf(dataset=list(pvalues=GWAS$P, N=GWAS$N, MAF=MAF, snp=GWAS$rsid,type="quant"))
  
  ##Perform coloc.abf test to test whether the causal variant is the same
  if (test_type=="P") {
    coloc <- coloc.abf(dataset1=list(pvalues=GWAS$P, N=GWAS$N, snp=GWAS$rsid, type="quant"),
                       dataset2=list(pvalues=sQTL$P, N=sQTL_N, snp=sQTL$rsid, type="quant"),
                       MAF=MAF)
  } else {
    coloc <- coloc.abf(dataset1=list(beta=GWAS$B, varbeta=GWAS$var, N=GWAS$N, snp=GWAS$rsid, type="quant"),
                       dataset2=list(beta=sQTL$B, varbeta=sQTL$var, N=sQTL_N, snp=sQTL$rsid, type="quant"),
                       MAF=MAF)
  }
  
  #Determine winning model. H4 means likeliest scenario is shared causal variant. H0=no causal variants, H1=only 1st dset has causal variant, 
  #H2=only 2nd dset has causal variant, H3=both have causal variants but they're different
  model<-which.max(coloc$summary[2:5])
  shared_var=FALSE
  if (model==4) {
    shared_var=TRUE
  }
  
  ##Plots
  sQTL_plot<-head(sQTL_finemap,-1)
  GWAS_plot<-head(GWAS_finemap,-1)
  
  sQTL_plot$genom_distance<-sQTL$genom_distance
  GWAS_plot$genom_distance<-GWAS$genom_distance
  
  #Get index of most probable causal SNP
  max_sQTL<-sQTL_plot[which.max(sQTL_plot$SNP.PP),]
  max_GWAS<-GWAS_plot[which.max(GWAS_plot$SNP.PP),]
  
  #Plot of -log10 p values of each
  Q<-ggplot(sQTL_plot)+
    geom_point(aes(x=genom_distance,y=SNP.PP),color="turquoise4")+
    ggtitle(paste(intron," sQTL p-values",sep=""))+
    ylab("log10 PP of SNP being causal")+
    xlab(paste("Chromosome ",chrom,sep=''))+
    annotate("point",x=max_sQTL$genom_distance,y=max_sQTL$SNP.PP,color="darkred",size=3)+
    annotate("text",x=max_sQTL$genom_distance,y=max_sQTL$SNP.PP,vjust=-0,hjust=-0.1,label=max_sQTL$snp,size=3)+  
    theme_minimal()
  
  G<-ggplot(GWAS_plot)+
    geom_point(aes(x=genom_distance,y=SNP.PP),color="steelblue4")+
    ggtitle(paste(trait," GWAS summary statistics",sep=""))+
    ylab("log10 PP of SNP being causal")+
    xlab(paste("Chromosome ",chrom,sep=''))+
    annotate("point",x=max_GWAS$genom_distance,y=max_GWAS$SNP.PP,color="darkred",size=3)+
    annotate("text",x=max_GWAS$genom_distance,y=max_GWAS$SNP.PP,vjust=-0,hjust=-0.1,label=max_GWAS$snp,size=3)+  
    theme_minimal()
  
  if (shared_var==TRUE) {
    #Plot of posterior probability of H4
    coloc_plot<-coloc$results
    coloc_plot$genom_distance<-sQTL$genom_distance
    max_coloc<-coloc_plot[which.max(coloc_plot$SNP.PP.H4),]
    
    C<-ggplot(coloc_plot)+
      geom_point(aes(x=genom_distance,y=SNP.PP.H4),color="darkgreen")+
      ggtitle(paste("Probability that SNP is causal for ",intron," sQTL and GWAS",sep=""))+
      ylab("log10 posterior probability")+
      xlab(paste("Chromosome ",chrom,sep=''))+
      annotate("point",x=max_coloc$genom_distance,y=max_coloc$SNP.PP.H4,color="darkred",size=3)+
      annotate("text",x=max_coloc$genom_distance,y=max_coloc$SNP.PP.H4,vjust=-0,hjust=-0.1,label=max_coloc$snp,size=3)+  
      annotate("text",x=quantile(coloc_plot$genom_distance,0.3),y=0.5,label=paste("PP of shared variant: ",coloc$summary[model],sep=''))+
      theme_minimal()
    
    pdf(paste("coloc/",intron,"_",trait,"_coloc_plots.pdf",sep=""),width=8,height=10)
    grid.arrange(Q,G,C,nrow=3)
  } else {
    pdf(paste("coloc/",intron,"_",trait,"_coloc_plots.pdf",sep=""),width=8,height=7)
    grid.arrange(Q,G,nrow=2)
  }
  
  
  write.table(sQTL_finemap,paste("coloc/sQTL_",intron,"_coloc_finemap",sep=""),quote=FALSE,row.names=FALSE,sep='\t')
  write.table(GWAS_finemap,paste("coloc/",intron,"_GWAS_coloc_finemap",sep=""),quote=FALSE,row.names=FALSE,sep='\t')
  write.table(coloc$summary,paste("coloc/",intron,"_coloc_summary",sep=""),quote=FALSE,row.names=TRUE,col.names=FALSE,sep='\t')
  write.table(coloc$results,paste("coloc/",intron,"_coloc_full",sep=""),quote=FALSE,row.names=FALSE,sep='\t')
  write("done",paste("coloc/",intron,"_coloc_status",sep=""))
} else {
  write("couln't run; empty locus files",paste("coloc/",intron,"_coloc_status",sep=""))
}




        
