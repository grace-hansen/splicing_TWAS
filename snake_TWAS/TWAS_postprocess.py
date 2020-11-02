import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("chrom", help="chromosome to analyze, in format 'chr1'")
parser.add_argument("tissue", help="tissue being analyzed")
parser.add_argument("data_source", help="origin of expression data (e.g. GTEx, CMC")
parser.add_argument("out", help="Output directory")
parser.add_argument("sumstats", help="GWAS summary stats")
parser.add_argument("LD", help="chromosome to analyze, in format 'chr1'")
args = parser.parse_args()

chromnum = args.chrom.split('r')[-1:][0]

##############TWAS postprocessing#################
#This scripts calculates significance for TWAS results, and does postprocessing on them (makes graphs, etc)
##################################################

#Calculate significance
numweights = sum(1 for f in os.listdir("%s/WEIGHTS"%(args.tissue)) if f.endswith(".wgt.RDat"))
cmd = "cat %s%s.%s.dat | awk 'NR == 1 || $NF < 0.05/%s' > %s%s.%s.dat.top" % (args.out, args.data_source, chromnum, numweights, args.out, args.data_source, chromnum)
print cmd
os.system(cmd)

#Make pretty graphs
cmd = """Rscript /project2/nobrega/grace/splicing/scripts/fusion_twas-master/FUSION.post_process.R \
--sumstats %s \
--tissue %s \
--input %s%s.%s.dat.top \
--out %s%s.%s.dat.top.analysis \
--ref_ld_chr %s \
--chr %s \
--plot \
--plot_individual \
--locus_win 100000""" % (args.sumstats, args.tissue, args.out, args.data_source, chromnum, args.out, args.data_source, chromnum, args.LD, chromnum)
print cmd
os.system(cmd)