#!/bin/bash

# prepares family for bcbio run when input files are family_sample.bam or family_sample_1/2.fq.gz
echo "Parameters:"
family=$1

echo "family="$family
# analysis:
# -default
# -noalign
# -align_decoy
# -cnvkit
# -sv: no align, bam has to be cleaned after aligning to decoy, call SV, use a bed file for SV regions
# -small_variants
# -validate
analysis=$2

echo "analysis="$analysis

cd $family

cp ~/cre/bcbio.sample_sheet_header.csv $family.csv

cd input

#there should be no other files except input fq.gz or bams in the input dir, bed files are ignored
ls | egrep -v ".bed$" | sed s/.bam// | sed s/.bai// | sed s/"_1.fq.gz"// | sed s/"_2.fq.gz"// | sort | uniq > ../samples.txt

cd ..

#no family is needed for single sample
while read sample
do
    echo $sample","$sample","$family",,," >> $family.csv
done < samples.txt

#default template
template=~/crg/config/crg.bcbio.default.yaml
if [ -n "$2" ]
then
    template=~/crg/config/crg.bcbio.$2.yaml
fi

echo "Using config template: " $template

bcbio_nextgen.py -w template $template $family.csv input/*

mv $family/config .
mv $family/work .
rm $family.csv
rmdir $family

cd ..
