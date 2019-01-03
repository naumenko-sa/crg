#!/bin/bash

#PBS -l walltime=00:30:00,nodes=1:ppn=2
#PBS -joe .
#PBS -d .
#PBS -l vmem=5g,mem=5g

TODAY=`date +%Y-%m-%d`
if [ -z $1 ]; then
	echo "Specify family ID as first arguement to script"
	exit
fi
FAMILY_ID=$1
OUT=${FAMILY_ID}.wgs.sv.${TODAY}.csv

HGMD=${HOME}/gene_data/HGMD_2018/hgmd_pro.db
EXON_BED=${HOME}/gene_data/protein_coding_genes.exons.fixed.sorted.bed
HPO=${HOME}/gene_data/HPO_2018/${FAMILY_ID}_HPO.txt
EXAC=${HOME}/gene_data/ExAC/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt
OMIM=${HOME}/gene_data/OMIM_2018-11-01/genemap2.txt

if [ ! $2 ]; then
	IN_FILES=`ls *.sv.csv | tr '\n' ' '`
else
	IN_FILES=${*:2}
fi

if [[ "$OSTYPE" == *"darwin"* ]]; then
	PY=python3
elif [[ "$OSTYPE" == "linux"* ]]; then
	module load python/3.5.6
	PY=python
fi

echo "${PY} ${HOME}/crg/crg.intersect_sv_reports.py -exon_bed=${EXON_BED} -hgmd=${HGMD} -hpo=${HPO} -exac=${EXAC} -omim=${OMIM} -o=${OUT} -i ${IN_FILES}"
${PY} ${HOME}/crg/crg.intersect_sv_reports.py -exon_bed=${EXON_BED} -hgmd=${HGMD} -hpo=${HPO} -exac=${EXAC} -omim=${OMIM} -o=${OUT} -i ${IN_FILES}

if [[ "$OSTYPE" == *"darwin"* ]]; then
	open -a 'Microsoft Excel' ${OUT}
fi
