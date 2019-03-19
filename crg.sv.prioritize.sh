#!/bin/bash
#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=21g,mem=21g

##################################################################################################
#   parameters:
#   $1 = case = project = family
#   $2 = panel = panel.bed - bed file with gene coordinates, created with ~/bioscripts/genes.R
#   $3 = dgv = tcag.tsv - report from TCAG which contains SV frequency in DGV_N_subjects columns
##################################################################################################

if [ -z $case ]
then
    case=$1
fi

if [ -z $panel ]
then
    panel=$2
fi

if [ -z $dgv ]
then
    dgv=$3
fi

#get only PASS calls
sample=`bcftools query -l ${case}-metasv.vcf.gz`

bcftools view -f PASS -o $sample.pass.vcf.gz -Oz ${case}-metasv.vcf.gz
tabix $sample.pass.vcf.gz

bedtools intersect -a $sample.pass.vcf.gz -b $panel -header -u > $sample.pass.region.vcf

#generating input file for TCAG annotation
crg.vcf2tsv.py $sample.pass.region.vcf > $sample.tsv

#svscore
if [ ! -f $sample.pass.region.svscore.vcf ]
then
    echo "Generating SV scores: " `date`
    SVSCORE_DATA=/hpf/largeprojects/ccmbio/arun/Tools/SVScore
    SVSCORE_SCRIPT=/hpf/largeprojects/ccmbio/naumenko/tools/svscore
    module load perl/5.20.1

    perl -w $SVSCORE_SCRIPT/svscore.pl -o max,sum,top5,top10,mean \
			 -e $SVSCORE_DATA/tests/refGene.exons.bed \
			 -f $SVSCORE_DATA/tests/refGene.introns.bed \
			 -dvc $SVSCORE_DATA/tests/whole_genome_SNVs.tsv.gz  \
			 -i $sample.pass.region.vcf > $sample.pass.region.svscore.vcf
fi

#generating final report with or without dgv
echo "Generating final report: " `date`
crg.sv.parse.py $sample.pass.region.svscore.vcf $dgv > $sample.sv.csv
