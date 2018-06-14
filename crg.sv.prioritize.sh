#!/bin/bash

###########################################################################################
#   parameters:
#   $1 = case = project = family
#   $2 = panel.bed - bed file with gene coordinates, created with ~/bioscripts/genes.R
#   $3 = tcag.tsv - report from TCAG which contains SV frequency in DGV_N_subjects columns
###########################################################################################

#get only PASS calls
sample=`bcftools query -l ${1}-metasv.vcf.gz`

bcftools view -f .,PASS -o $sample.pass.vcf.gz -Oz ${1}-metasv.vcf.gz
tabix $sample.pass.vcf.gz

bedtools intersect -a $sample.pass.vcf.gz -b $2 -header -u > $sample.pass.region.vcf

crg.sv.parse.py $sample.pass.region.vcf $3 > $sample.sv.csv
