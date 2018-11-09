#!/bin/bash

#PBS -l walltime=00:30:00,nodes=1:ppn=2
#PBS -joe .
#PBS -d .
#PBS -l vmem=5g,mem=5g

TODAY=`date +%Y-%m-%d`
FAMILY_ID=$1
OUT=${FAMILY_ID}.wgs.sv.${TODAY}.csv

if [[ "$OSTYPE" == "linux"* ]]; then
	HGMD=${HOME}/gene_data/HGMD_2018/hgmd_pro.db
	EXON_BED=${HOME}/gene_data/protein_coding_genes.exons.fixed.sorted.bed
fi

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

echo "${PY} ${HOME}/crg/crg.intersect_sv_reports.py -exon_bed=${EXON_BED} -hgmd=${HGMD} -o=${OUT} -i ${IN_FILES}"
${PY} ${HOME}/crg/crg.intersect_sv_reports.py -exon_bed=${EXON_BED} -hgmd=${HGMD} -o=${OUT} -i ${IN_FILES}

if [[ "$OSTYPE" == *"darwin"* ]]; then
	open -a 'Microsoft Excel' ${OUT}
fi
