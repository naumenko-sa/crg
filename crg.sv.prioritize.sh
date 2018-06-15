#!/bin/bash

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=21g,mem=21g

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

#svscore

if [ ! -f $sample.pass.region.svscore.vcf ]
then

    SVSCORE_DATA=/hpf/largeprojects/ccmbio/arun/Tools/SVScore
    SVSCORE_SCRIPT=/hpf/largeprojects/ccmbio/naumenko/tools/svscore
    module load perl/5.20.1

    perl -w $SVSCORE_SCRIPT/svscore.pl -o max,sum,top5,top10,mean \
			 -e $SVSCORE_DATA/tests/refGene.exons.bed \
			 -f $SVSCORE_DATA/tests/refGene.introns.bed \
			 -dvc $SVSCORE_DATA/tests/whole_genome_SNVs.tsv.gz  \
			 -i $sample.pass.region.vcf > $sample.pass.region.svscore.vcf
fi

crg.sv.parse.py $sample.pass.region.svscore.vcf $3 > $sample.sv.csv
