#!/bin/bash

#####################
# parameters:
#   $1 = case = project = family
#   $2 = panel.bed - bed file with gene coordinates, created with ~/bioscripts/genes.R

#get only PASS calls
sample=`bcftools query -l ${1}-metasv.vcf.gz`

bcftools view -f .,PASS -o $sample.pass.vcf.gz -Oz ${1}-metasv.vcf.gz
tabix $sample.pass.vcf.gz

bedtools intersect -a $sample.pass.vcf.gz -b $2 -header > $sample.pass.region.vcf
