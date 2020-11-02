import os, glob, argparse, commands, gzip, shutil, mmap, multiprocessing as mp
import pandas as pd
parser = argparse.ArgumentParser()
parser.add_argument("tissue", help="tissue being analyzed")
parser.add_argument("trait", help="tissue being analyzed")
parser.add_argument("chrom", help="chromosome to analyze, in format 'chr1'")
parser.add_argument("data_source", help="origin of expression data (e.g. GTEx, CMC")
parser.add_argument("offset", help="size of offset in scientific notation, e.g. 5e4")
parser.add_argument("covars", help="Covariates to include when making TWAS weights")
parser.add_argument("temp", help="Location of temp directory for snakefile flag")
args = parser.parse_args()

##############Make plinks and compute_TWAS_weights#################
# This  script takes genotype and phenotype(here, intron counts in *perind.counts format) data from GTEx and prepares plinks for TWAS.
# It then takes these plinks and  computes TWAS weights.
###################################################################

######### Chunking function for breaking the file up into parts for multiprocessing ######
def chunk(max, num):
    seq=range(max)
    avg = len(seq) / float(num)
    out = []
    last = 0.0
    while last < len(seq):
        out.append([int(last),int(last + avg)])
        last += avg
    return out
###########
chromnum = args.chrom.split('r')[-1:][0]
tissue=args.tissue
trait=args.trait
offset = float(args.offset)
temp=args.temp
data_source=args.data_source
covars=args.covars

#Make directories to use later
tmp = "%s/tmp_%s" % (tissue,chromnum)
if not os.path.isdir(tmp):
	os.mkdir(tmp)

weightsdir = "%s/WEIGHTS_%s" % (tissue,chromnum)
if not os.path.isdir(weightsdir):
	os.mkdir(weightsdir)

#Make list of all variants in GWAS summary statistics to filter by
if not os.path.exists("GWAS_rsids_%s"%(chromnum)):
	cmd="cut -f2 /project2/nobrega/grace/genos/%s/plink/chr%s_rsids.bim"%(data_source,chromnum)
	status,rsids=commands.getstatusoutput(cmd)
	rsids=list(rsids.split('\n'))
	GWAS=pd.read_csv("%s_GWAS_sumstats"%(trait),sep='\t')
	with open("GWAS_rsids_%s"%(chromnum),'w') as file:
		GWAS_rsids=GWAS["SNP"].tolist()
		intersect=list(set(rsids).intersection(GWAS_rsids))
		intersect='\n'.join(intersect)
		file.write(intersect+'\n')
	del GWAS

GWAS_rsid_file="GWAS_rsids_%s"%(chromnum)

subject_IDs="%s/%s_subject_IDs"%(tissue,tissue)

######## Start making weights
#Break intron file into chunks
num_threads=16
introns=gzip.open('%s/%s_perind.decs.%s.gz'%(tissue,tissue,chromnum))
num_lines = sum(1 for line in introns)
introns.close()
chunks=chunk(num_lines,num_threads)

introns=gzip.open('%s/%s_perind.decs.%s.gz'%(tissue,tissue,chromnum),'r')
intron_mat=introns.readlines()
intron_mat=intron_mat[1:] 
introns.close()

#For some reason, multiprocessing DOES NOT LIKE MULTIPLE OPEN COMPRESSED FILES and throws lots of very strange errors. Therefore, we read all the lines in, and then close the file.

######### Open file record for whether chunk has been run
if not os.path.exists("%s/clusts_done_%s.txt"%(tissue,chromnum)):
	record=open("%s/clusts_done_%s.txt"%(tissue,chromnum),'w+')
	s=''
elif os.stat("%s/clusts_done_%s.txt"%(tissue,chromnum)).st_size == 0:
	record=open("%s/clusts_done_%s.txt"%(tissue,chromnum),'w+')
	s=''
else:
	record=open("%s/clusts_done_%s.txt"%(tissue,chromnum),'a+')
	s = mmap.mmap(record.fileno(), 0, access=mmap.ACCESS_READ)
	record.close()

def make_weights(linerange):
	i=0
	chunk_record=open("%s/clusts_done_%s_%s.txt"%(tissue,chromnum,linerange[0]),'w')
	for ln in intron_mat:
		if i >= linerange[0] and i <= linerange[1]:
			clust=ln.split()[:1]
			clust[0]=clust[0].decode("utf-8")
			if s.find(clust[0]) == -1:
				chunk_record.write(clust[0]+'\n')
				c,start,end,clust_ID=clust[0].split(":")[0:4]
				if c[0]=='c':
					c=c.split('r')[1]
				clust=clust_ID+'-'+c+'-'+start+'-'+end
				if data_source=="CMC"  or data_source=="GTEx_v8":
					vals = ln.split()[1:]
					vals=[v.decode("utf-8") for v in vals]
				if data_source=="GTEx_v7":
					vals = ln.split()[3:]
				#Make phenotype file of subject IDs + intron prevalence to filter plink by
				with open("%s_clust"%(clust),'w') as file:
					file.write('\n'.join(vals)+'\n')
				clust_vals="%s_clust"%(clust)
				cmd="paste %s %s > %s_cluster_prevalence"%(subject_IDs,clust_vals,clust)
				os.system(cmd)
				clust_vals_file="%s_cluster_prevalence"%(clust)
				#Define window within which to grab variants
				p0=int(start)-offset
				if p0 < 0:
					p0 = 0
				p1=int(end)+offset
				#Create plinks to run compute_weights on
				cmd="plink --bfile /project2/nobrega/grace/genos/%s/plink/chr%s_rsids \
				--pheno %s --make-bed --out %s/%s --keep %s \
				--chr %s --from-bp %s --to-bp %s --extract %s"%(data_source,chromnum,clust_vals_file,tissue,clust,clust_vals_file,chromnum,p0,p1,GWAS_rsid_file)
				os.system(cmd)
				os.remove("%s_clust"%(clust))
				os.remove("%s_cluster_prevalence"%(clust))
				#Run compute_weights, delete plink and temp files after
				if c == chromnum:
					if covars=="None":
						cmd = """Rscript /project2/nobrega/grace/splicing/scripts/fusion_twas-master/FUSION.compute_weights.R \
						--bfile %s/%s \
						--tmp %s/tmp_%s/ \
						--save_hsq \
						--out %s/WEIGHTS_%s/%s \
						--PATH_gcta /project2/nobrega/grace/splicing/scripts/gcta_1.91.4beta/gcta64 \
						--models top1,blup,lasso,enet""" % (tissue, clust, tissue, chromnum, tissue, chromnum, clust)
					else:
						cmd = """Rscript /project2/nobrega/grace/expression/scripts/fusion_twas-master/FUSION.compute_weights.R \
						--bfile %s/%s \
						--covar %s \
						--tmp %s/tmp_%s/ \
						--save_hsq \
						--out %s/WEIGHTS_%s/%s \
						--PATH_gcta /project2/nobrega/grace/expression/scripts/gcta_1.91.4beta/gcta64 \
						--models top1,blup,lasso,enet""" % (tissue, clust, covars, tissue, chromnum, tissue, chromnum, clust)
					os.system(cmd)
					os.system("cat %s/WEIGHTS_%s/%s.hsq >> %s/WEIGHTS_%s/chr_%s.hsq && rm %s/WEIGHTS_%s/%s.hsq"%(tissue, chromnum, clust, tissue, chromnum, chromnum, tissue, chromnum, clust))
				for p in glob.glob("%s/%s.*"%(tissue,clust)):
					os.remove(p)
		if i > linerange[1]:
			chunk_record.close()
			cmd="cat %s/clusts_done_%s_%s.txt >> %s/clusts_done_%s.txt && rm %s/clusts_done_%s_%s.txt"%(tissue,chromnum,linerange[0],tissue,chromnum,tissue,chromnum,linerange[0])
			os.system(cmd)
			break
		i=i+1
	chunk_record.close()
	cmd="cat %s/clusts_done_%s_%s.txt >> %s/clusts_done_%s.txt && rm %s/clusts_done_%s_%s.txt"%(tissue,chromnum,linerange[0],tissue,chromnum,tissue,chromnum,linerange[0])
	os.system(cmd)

pool=mp.Pool(processes=num_threads)
results=[pool.apply(make_weights,args=(chunk,)) for chunk in chunks]

#Remove temp directories if empty
shutil.rmtree(tmp,ignore_errors=True)

#Create temp file indicating TWAS weights have been made
temp="%smake_weights_%s.txt" % (temp,chromnum)
open(temp,'a').close()