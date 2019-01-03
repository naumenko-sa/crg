#!/bin/bash

# https://academic.oup.com/bioinformatics/article/33/7/1083/2748212
# https://github.com/lganel/SVScore

echo "Generating SV scores: " `date`
SVSCORE_DATA=/hpf/largeprojects/ccmbio/arun/Tools/SVScore
SVSCORE_SCRIPT=/hpf/largeprojects/ccmbio/naumenko/tools/svscore
module load perl/5.20.1

sample=$1

perl -w $SVSCORE_SCRIPT/svscore.pl -o max,sum,top5,top10,mean \
		 -e $SVSCORE_DATA/tests/refGene.exons.bed \
		 -f $SVSCORE_DATA/tests/refGene.introns.bed \
		 -dvc $SVSCORE_DATA/tests/whole_genome_SNVs.tsv.gz  \
		 -i $sample.pass.region.vcf > $sample.pass.region.svscore.vcf
