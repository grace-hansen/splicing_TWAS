import os, subprocess, glob, argparse
parser = argparse.ArgumentParser()
parser.add_argument("chrom", help="chromosome to analyze, in format 'chr1'")
parser.add_argument("tissue", help="tissue being analyzed")
parser.add_argument("data_source", help="origin of expression data (e.g. GTEx, CMC")
parser.add_argument("out", help="Output directory")
parser.add_argument("sumstats", help="GWAS summary stats")
parser.add_argument("LD", help="chromosome to analyze, in format 'chr1'")
args = parser.parse_args()

##############Run TWAS#################
#This script takes weights in a chromosome-specific WEIGHTS_* directory and runs TWAS on them.
#######################################

chromnum = args.chrom.split('r')[-1:][0]

#Move all chromosome-specific WEIGHTS_* dir contents to WEIGHTS
if not os.path.isdir("%s/WEIGHTS"%(args.tissue)):
	os.mkdir("%s/WEIGHTS"%(args.tissue))
cmd = "mv %s/WEIGHTS_%s/* %s/WEIGHTS && rm -r %s/WEIGHTS_%s" % (args.tissue,chromnum,args.tissue,args.tissue,chromnum)
os.system(cmd)

#make WEIGHTS*.pos file listing all *RDat files and corresponding gene coordinates
posfile="%s/WEIGHTS/%s.%s.pos" % (args.tissue,args.data_source,chromnum)
pos = open(posfile,'w')
pos.write("WGT\tID\tCHR\tP0\tP1\n")

weights=glob.glob("%s/WEIGHTS/*RDat"%(args.tissue))
for ln in weights:
	line = ln.strip().split('/')[2]
	clust = ln.strip().split('.wgt.RDat')[0].split('/')[2]
	clustID,c,s,e = clust.split('-')
	if c == chromnum:
		posline=("%s\t%s\t%s\t%s\t%s\n") % (line, clust, chromnum, s, e)
		pos.write(posline)
pos.close()

#Run TWAS
cmd = """Rscript /project2/nobrega/grace/splicing/scripts/fusion_twas-master/FUSION.assoc_test.R \
--sumstats %s \
--weights ./%s/WEIGHTS/%s.%s.pos \
--weights_dir ./%s/WEIGHTS/ \
--ref_ld_chr %s \
--chr %s \
--out %s%s.%s.dat""" % (args.sumstats, args.tissue, args.data_source, chromnum, args.tissue, args.LD, chromnum, args.out, args.data_source, chromnum)
os.system(cmd)