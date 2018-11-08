#!/bin/bash

#PBS -l walltime=00:30:00,nodes=1:ppn=2
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

TODAY=`date +%Y-%m-%d`
FAMILY_ID=$1
OUT=${FAMILY_ID}.wgs.sv.${TODAY}.csv
HGMD=/Users/denniskao/hgmd/hgmd_pro.db

if [ ! $2 ]; then
	IN_FILES=`ls *.sv.csv | tr '\n' ' '`
else
	IN_FILES=${*:2}
fi

if [[ "$OSTYPE" == *"darwin"* ]]; then
	PY=python3
elif [[ "$OSTYPE" == "linux"* ]]; then
	module load python/3.5.2
	PY=python
fi

echo "${PY} /Users/denniskao/crg/crg.intersect_sv_reports.py -exon_bed=/Users/denniskao/protein_coding_genes.exons.fixed.sorted.bed -hgmd=${HGMD} -o=${OUT} -i ${IN_FILES}"
${PY} /Users/denniskao/crg/crg.intersect_sv_reports.py -exon_bed=/Users/denniskao/protein_coding_genes.exons.fixed.sorted.bed -hgmd=${HGMD} -o=${OUT} -i ${IN_FILES}

if [[ "$OSTYPE" == *"darwin"* ]]; then
	open -a 'Microsoft Excel' ${OUT}
fi