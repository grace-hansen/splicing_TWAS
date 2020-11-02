import os
import argparse
import commands
import gzip
import pandas as pd
parser = argparse.ArgumentParser()
parser.add_argument("tissue", help="tissue being analyzed")
parser.add_argument("trait", help="trait being analyzed")
parser.add_argument("chrom", help="chromosome to analyze, in format 'chr1'")
parser.add_argument("data_source", help="data source, e.g. CMC, GTEx")
parser.add_argument("temp", help="Location of temp directory for snakefile flag")
parser.add_argument("vcf", help="Location of genotype vcf")
parser.add_argument("--exclude", help="List of subject IDs to exclude (e.g. non-Europeans",default="None")
args = parser.parse_args()

##############Keep only IDs with genotype info#################
chromLoc = args.chrom
chromnum = chromLoc.split('r')[-1:][0]
trait=args.trait
tissue=args.tissue
temp=args.temp
data_source=args.data_source
vcf=args.vcf

#Convert perind.counts into pandas
raw_counts=pd.read_csv("/project2/nobrega/grace/splicing/%s/%s/%s_perind.counts.gz"%(trait,tissue,tissue),sep=' ',compression='gzip')

#Get IDs from perind.counts file
splice_IDs=list(raw_counts.columns.values)

#Get IDs from genotype files
if data_source == "GTEx_v8" or data_source == "GTEx_v7":
	cmd="bcftools query -l %s"%(vcf)
	status, output = commands.getstatusoutput(cmd)
	geno_IDs=list(output.split('\n'))
if data_source =="CMC":
	cmd="cut -f1 -d' ' /project2/nobrega/grace/genos/CMC/plink/chr1_rsids.fam"
	status, output = commands.getstatusoutput(cmd)
	geno_IDs=list(output.split('\n'))

#Filter perind.counts IDs, keeping only those in genotype.vcf
new_splice_IDs=list(set(splice_IDs) & set(geno_IDs))
print new_splice_IDs

#Filter out excluded IDs, if given
if args.exclude != "None":
    with open(args.exclude,'r') as excl_file:
        exclude_IDs=excl_file.readlines()
        exclude_IDs = ' '.join(exclude_IDs).replace('\n','').split()
    new_splice_IDs=list(set(new_splice_IDs)-set(exclude_IDs))

new_splice_IDs= ["chrom"] + new_splice_IDs
counts=raw_counts[new_splice_IDs]


#Write out list of subject IDs if it hasn't been written yet
if not os.path.exists("/project2/nobrega/grace/splicing/%s/%s/%s_subject_IDs"%(trait,tissue,tissue)):
	subj_list = open("/project2/nobrega/grace/splicing/%s/%s/%s_subject_IDs"%(trait,tissue,tissue), "w")
	for ID in new_splice_IDs:
		if ID != "chrom":
			subj_list.write(ID+'\t'+ID+'\n')
	subj_list.close()


##############Convert ratios to decimal values#################
out=gzip.open("/project2/nobrega/grace/splicing/%s/%s/%s_perind.decs.%s.gz"%(trait,tissue,tissue,chromnum),'w')
IDs=list(counts.columns.values)
IDs.remove(IDs[0])
out.write('clust'+'\t'+'\t'.join(IDs)+'\n')

for ln in counts.itertuples():
	header=ln[1]
	c,s,e,clust_ID=header.split(":")[0:4]
	if c[0]=='c':
		c=c.split('r')[1]
	if c == str(chromnum):
		row=header
		vals=ln[2:]
		for i in vals:
			num,den=i.split('/')
			try:
				dec=round(float(num)/float(den),4)
			except ZeroDivisionError:
				dec=0
			row=row+'\t'+str(dec)
		row=row+'\n'
		out.write(row)
out.close()