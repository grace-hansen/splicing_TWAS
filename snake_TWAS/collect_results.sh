#! /bin/sh

trait=$1
tissue=$2
data_source=$3

cd /project2/nobrega/grace/splicing/$trait/$tissue/results

head -1 ${data_source}.1.dat.top > ${data_source}.all.dat.top
for i in $(seq 1 22); do
    sed '1d' ${data_source}.$i.dat.top >> temp
done
sort -g -k19,19 temp >> ${data_source}.all.dat.top
rm temp

head -1 ${data_source}.1.dat > ${data_source}.all.dat
for i in $(seq 1 22); do
    sed '1d' ${data_source}.$i.dat >> temp
done
sort -g -k19,19 temp >> ${data_source}.all.dat
rm temp

if [ ${data_source} == "GTEx_v7" ]; then
	python2 /project2/nobrega/grace/splicing/scripts/get_results_gene.py /project2/nobrega/grace/gencode.v19.annotation_hg19.gtf ${data_source}.all.dat.top
	python2 /project2/nobrega/grace/splicing/scripts/get_results_gene.py /project2/nobrega/grace/gencode.v19.annotation_hg19.gtf ${data_source}.all.dat
elif  [ ${data_source} == "CMC" ]; then
	python2 /project2/nobrega/grace/splicing/scripts/get_results_gene.py /project2/nobrega/grace/gencode.v19.annotation_hg19.gtf ${data_source}.all.dat.top
	python2 /project2/nobrega/grace/splicing/scripts/get_results_gene.py /project2/nobrega/grace/gencode.v19.annotation_hg19.gtf ${data_source}.all.dat
else 
	python2 /project2/nobrega/grace/splicing/scripts/get_results_gene.py /project2/nobrega/grace/gencode.v29.annotation_hg38.gtf ${data_source}.all.dat.top
	python2 /project2/nobrega/grace/splicing/scripts/get_results_gene.py /project2/nobrega/grace/gencode.v29.annotation_hg38.gtf ${data_source}.all.dat
fi

mkdir posthoc
sed '1d' ${data_source}.all.dat.top.gene | cut -f1,20 | sort -u -k1,1 | sort -g -k2,2 | cut -f1 > posthoc/sig_genes
sed '/^NA$/d' posthoc/sig_genes > temp
mv temp posthoc/sig_genes

sed '1d' ${data_source}.all.dat.top.gene | cut -f2 | cut -f5  -d'/' | cut -f1 -d'.' > introns
sed '1d' ${data_source}.all.dat.top.gene | cut -f20 > pvals
paste introns pvals | sort -u -k1,1 | sort -g -k2,2 | cut -f1 > posthoc/sig_introns
rm introns pvals