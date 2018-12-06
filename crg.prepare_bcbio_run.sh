#!/bin/bash

# prepares family for bcbio run when input files are family_sample.bam or family_sample_1/2.fq.gz
family=$1

# template type = default | noalign | align_decoy | cnvkit
noalign=$2

cd $family

cp ~/cre/bcbio.sample_sheet_header.csv $family.csv

cd input

#there should be no other files except input fq.gz or bams in the input dir
ls | sed s/.bam// | sed s/.bai// | sed s/"_1.fq.gz"// | sed s/"_2.fq.gz"// | sort | uniq > ../samples.txt

cd ..

#no family is needed for single sample
while read sample
do
    echo $sample","$sample","$family",,," >> $family.csv
done < samples.txt

#default template
template=~/crg/crg.templates.default.yaml

template_type=$2

if [ -n "$2" ]
then
    if [ $template_type == "noalign" ]
    then
	template=~/crg/crg.templates.noalign.yaml
    elif [ $template_type == "align_decoy" ]
    then
	template=~/crg/crg.templates.align_decoy.yaml
    else
	template=~/crg/crg.templates.default.yaml
    fi
fi

bcbio_nextgen.py -w template $template $family.csv input/*

mv $family/config .
mv $family/work .
rm $family.csv
rmdir $family

cd ..
