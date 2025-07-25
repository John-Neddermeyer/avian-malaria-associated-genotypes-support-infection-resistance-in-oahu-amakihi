#!/bin/bash

REF=/projects/foster_lab/john/birds/oaam/oaam.v1.simple.fa
conda activate bcftools_env
srun bcftools mpileup --threads 24 --max-depth 100 -f $REF \
-Q 20 -q 20 --bam-list list_of_association_bam_files.txt --output-type v | \
bcftools call --threads 24 -m -O z -f GQ \
-o association.bcftools.called.vcf.gz

bcftools filter -g3 -G10 -e'QUAL<20 \
|| (RPBZ<0.1 && QUAL<25) || (AC<2 && QUAL<25) || MAX(IDV)<=3 || \
MAX(IDV)/MAX(DP)<=0.3' association.bcftools.called.vcf.gz -O z \
-o filtered_bcftools_called.vcf.gz --threads 12 --write-index

VCFTOOLS=/scratch/jhn43/tools/vcftools/bin/vcftools
VCF_FILE=filtered_bcftools_called.vcf.gz
OUT_VCF=oaam_full_filtered.vcf.gz

$VCFTOOLS --gzvcf $VCF_FILE --remove-indels --maf 0.10 --max-missing 1 \
--min-alleles 2 --max-alleles 2 --recode --stdout | gzip > $OUT_VCF
