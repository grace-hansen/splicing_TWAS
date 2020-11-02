#!/bin/bash
#This script makes an LD reference panel from GTEx v7 data using LDAK version 5..

LD_dir=$1
chrom=$2
data_source=$3
vcf=$4
lookup=$5
excludes=$6
MAF=$7


mkdir -p $LD_dir
cd $LD_dir

if [[ $data_source == *"GTEx"* ]]; then
	#Replace GTEx varIDs with rsids (only relevant for TWAS LD panel)
	if [ ! -e ~/midway/genos/${data_source}/vcf/chr${chrom}_rsid.vcf.gz ]; then
		python2 ~/midway/expression/scripts/add_rsids_raw_vcfs.py chr${chrom} ${data_source} $vcf $lookup
	fi

	#Keep only europeans (list made from grab_race_IDs.R in scripts/), subset by MAF
	bcftools query -l $vcf | tr ' ' '\n' > IDs_${chrom}
	comm -13 <(sort $excludes) <(sort IDs_${chrom}) > tmp_${chrom}
	paste tmp_${chrom} tmp_${chrom} > IDs.keep_${chrom}
	plink --vcf ~/midway/genos/${data_source}/vcf/chr${chrom}_rsid.vcf.gz --make-bed --out chr${chrom} --maf $MAF --keep IDs.keep_${chrom}
	echo "chr${chrom}" > chr${chrom}.txt

else
	cut -f1 -d' ' ~/midway/genos/${data_source}/plink/chr${chrom}_rsids.fam > IDs_${chrom}
	comm -13 <(sort $excludes) <(sort IDs_${chrom}) > tmp_${chrom}
	paste tmp_${chrom} tmp_${chrom} > IDs.keep_${chrom}
	plink --bfile ~/midway/genos/${data_source}/plink/chr${chrom}_rsids --make-bed --out chr${chrom} --maf $MAF --keep IDs.keep_${chrom}
	echo "chr${chrom}" > chr${chrom}.txt
fi

#Incorporate heritability info with LDAK
/project2/nobrega/grace/expression/scripts/ldak5.linux --make-bed chr${chrom}_LD --mbfile chr${chrom}.txt --exclude-odd YES --exclude-dups YES

#Incorporate genetic distances
plink --bfile chr${chrom}_LD --chr ${chrom} --cm-map  ~/midway/genos/genetic_map_b37/genetic_map_chr@_combined_b37.txt --make-bed --out map${chrom}
cat map${chrom}.bim | awk '{print $2, $3}' > map${chrom}.tmp

#Put into output format
awk '(NR==FNR){arr[$1]=$2;next}{print $1, $2, arr[$2], $4, $5, $6}' map${chrom}.tmp chr${chrom}_LD.bim > ${chrom}.bim
cp chr${chrom}_LD.fam ${chrom}.fam
cp chr${chrom}_LD.bed ${chrom}.bed

rm tmp_${chrom}
rm IDs_${chrom}
rm IDs.keep_${chrom}
rm chr${chrom}*
rm map${chrom}*
