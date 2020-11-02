import gzip, os, argparse, subprocess
import pandas as pd
parser = argparse.ArgumentParser()
parser.add_argument("trait", help="trait (name of directory in which tissue is located")
parser.add_argument("tissue", help="tissue")
parser.add_argument("intron", help="intron for which to extract information. Must be formatted like clu_41766-2-25174423-25176419")
parser.add_argument("GWAS_sumstats", help="full path to sumstats that will be used for coloc. Will only include SNPs in GWAS sumstats. rsid must be first column")
parser.add_argument("LD", help="LD reference panel being used (for ref/alt allele matching, in the format of 'LD")
parser.add_argument("data_source", help="data_source,e.g. CMC, GTEx")


##############Grab intron rsids#################
args = parser.parse_args()
trait=args.trait
tissue=args.tissue
intron=args.intron
GWAS_sumstats=args.GWAS_sumstats
LD=args.LD
data_source=args.data_source

chrom=intron.split('-')[1]

def make_intron_file(intron,trait,tissue):
    #Given the fastqtl file for all tested sQTL-intron relationships in a specific tissue, extract only the ones for your intron of interest,
        # associate the intron with a gene and include that in the output.
    clust_id,clust_chrom,intron_start,intron_stop=intron.split("-")
    print("writing sQTL entries for %s..."%(intron))
    with gzip.open(
        "/project2/nobrega/grace/splicing/%s/%s/%s_sQTLs_%s_%s.txt.gz"%(trait,tissue,tissue,intron,data_source),'wt') as out, gzip.open(
        "/project2/nobrega/grace/splicing/%s/%s/%s_sQTLs_%s.txt.gz"%(trait,tissue,tissue,data_source),'rt') as sQTL_dat:
        if data_source=="CMC":
            sQTL_dat.readline()
        for ln in sQTL_dat: #Put intron's vcf here
            if data_source == "GTEx_v7":
                sQTL_intron,var_id,dist,P,B= ln.split()[:5]
                sQTL_chrom,sQTL_start,sQTL_stop,sQTL_clust_id=sQTL_intron.split(":")[0:4]
                if sQTL_clust_id==clust_id and sQTL_start==intron_start and sQTL_stop==intron_stop:
                    out.write(ln)
            if data_source == "CMC":
                sQTL_intron,var_id,dist,P,B= ln.split()[:5]
                sQTL_chrom,sQTL_start,sQTL_stop,sQTL_clust_id=sQTL_intron.split(":")[0:4]
                if sQTL_clust_id==clust_id:
                    if sQTL_start==intron_start and sQTL_stop==intron_stop:
                        out.write(ln)


def merge_GWAS_sQTL(GWAS_sumstats,intron,trait,tissue):
    #Inner join the file made by make_intron_file and a supplied GWAS summary statistics file, so that the same rsids are in both. 
        #Save out these files
    #Load GWAS sumstats and sQTL sumstats
    GWAS=pd.read_csv(GWAS_sumstats,sep='\t')
    GWAS.rename(columns={ GWAS.columns[0]: "rsid" },inplace=True)
    GWAS_cols=list(GWAS.columns)
    try:
        sQTL=pd.read_csv("/project2/nobrega/grace/splicing/%s/%s/%s_sQTLs_%s_%s.txt.gz"%(trait,tissue,tissue,intron,data_source),delim_whitespace=True,compression='gzip',header=None)
    except: #No columns in sQTL = no variants within range of intron to test
        with open("/project2/nobrega/grace/splicing/%s/%s/results/posthoc/coloc/loci/%s_locus_status"%(trait,tissue,intron),'w') as status:
            status.write("not run; no sQTLs within locus")
        os.remove("/project2/nobrega/grace/splicing/%s/%s/%s_sQTLs_%s_%s.txt.gz"%(trait,tissue,tissue,intron,data_source))
        return
    if data_source=="GTEx_v7":
        sQTL.columns=["intron_id","variant_id","distance","P","B"]
        #Load up GTEx lookup, add reference allele and rsid columns to sQTL data
        sQTL_ids=sQTL['variant_id'].tolist()
        lookup=pd.read_csv("/project2/nobrega/GTEx/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz",sep='\t',compression="gzip")
        sQTL_lookup=lookup[lookup['variant_id'].isin(sQTL_ids)]
        sQTL_lookup.rename(columns={ "rs_id_dbSNP147_GRCh37p13": "rsid" },inplace=True)
        sQTL=sQTL.merge(sQTL_lookup[['rsid','variant_id','ref']],on='variant_id')
        del lookup
        #Load up GTEx eQTL data, extract MAF info from it
        MAF_dat=pd.read_csv("~/midway/GTEx_vars_MAF_chr%s.txt"%(chrom),sep='\t',header=None)
        MAF_dat.columns=['chr','loc','variant_id','MAF_ref','MAF']
        MAF_dat=MAF_dat[MAF_dat['variant_id'].isin(sQTL_ids)]
        sQTL=sQTL.merge(MAF_dat[['variant_id','MAF_ref','MAF']],on='variant_id')
        del MAF_dat
    if data_source=="CMC":
        sQTL.columns=["intron_id","rsid","distance","P","B","ref","MAF"]
        sQTL["MAF_ref"]=sQTL["ref"]
    merged=pd.merge(GWAS,sQTL,on='rsid')
    print merged.columns
    sQTL_cols=['intron_id', 'rsid','distance', 'P', 'B', 'ref','MAF'] #Want MAF and rsid in here, but not variant_id
    ##Flip alleles and MAF
    print "Flipping alleles so eQTL, GWAS, and LD panel match..."
    flips=[]
    for row in merged.itertuples():
        new_row=list(row)
        rsid=row.rsid
        #Get LD reference allele
        cmd="""awk -F ' ' '{ if ($2 == "%s") print $6}' %s*.bim"""%(rsid,LD)
        LD_ref=subprocess.check_output(cmd,shell=True).decode("utf-8").rstrip()
        if LD_ref != "":
        #Get sQTL reference allele
            sQTL_ref=row.ref
            if sQTL_ref != LD_ref:
                B=float(new_row[10])
                new_row[10]=-B
        #Get GWAS reference allele
            cmd="""awk -F '\t' '{ if ($1 == "%s") print $3}' %s """%(rsid,GWAS_sumstats)
            GWAS_ref=subprocess.check_output(cmd,shell=True).decode("utf-8").rstrip()
            if GWAS_ref != LD_ref:
                A1, A2, Z = new_row[2],new_row[3],float(new_row[4]) #Flip A1 and A2
                new_row[2]=A2
                new_row[3]=A1
                new_row[4]=-Z
        #Flip MAF
            MAF_ref=row.MAF_ref
            if MAF_ref != LD_ref:
                new_row[13]=1-float(row.MAF)
            new_row=new_row[1:]
            flips.append(new_row)
    merged_flipped=pd.DataFrame(flips,columns=merged.columns)
    #Keep only GWAS sumstats where rsid is in sQTL dataset
    GWAS_out=merged_flipped.ix[:,GWAS_cols]
    GWAS_out.to_csv("/project2/nobrega/grace/splicing/%s/%s/results/posthoc/coloc/loci/%s_GWAS_sumstats"%(trait,tissue,intron),sep='\t',index=False)
    #Keep only sQTL sumstats where rsid is in GWAS dataset
    sQTL_out=merged_flipped.ix[:,sQTL_cols]
    sQTL_out.to_csv("/project2/nobrega/grace/splicing/%s/%s/results/posthoc/coloc/loci/%s_sQTL.txt"%(trait,tissue,intron),sep='\t',index=False)
    with open("/project2/nobrega/grace/splicing/%s/%s/results/posthoc/coloc/loci/%s_locus_status"%(trait,tissue,intron),'w') as status:
            status.write("done")
    os.remove("/project2/nobrega/grace/splicing/%s/%s/%s_sQTLs_%s_%s.txt.gz"%(trait,tissue,tissue,intron,data_source))


try:
    os.makedirs("coloc/loci/")
except:
    pass


if not os.path.isfile("/project2/nobrega/grace/splicing/%s/%s/results/posthoc/coloc/loci/%s_sQTL.txt"%(
    trait,tissue,intron)) and not os.path.isfile("/project2/nobrega/grace/splicing/%s/%s/results/posthoc/coloc/loci/%s_GWAS_sumstats"%(
    trait,tissue,intron)):
    if not os.path.isfile("/project2/nobrega/grace/splicing/%s/%s/%s_sQTL_%s_rsid.txt.gz"%(trait,tissue,tissue,intron)):
        make_intron_file(intron,trait,tissue)
        merge_GWAS_sQTL(GWAS_sumstats,intron,trait,tissue)
    else: 
        merge_GWAS_sQTL(GWAS_sumstats,intron,trait,tissue)
