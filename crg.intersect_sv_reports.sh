#!/bin/bash

#PBS -l walltime=00:30:00,nodes=1:ppn=2
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

TODAY=`date +%Y-%m-%d`
FAMILY_ID=$1
OUT=${FAMILY_ID}.wgs.sv.${TODAY}.csv

if [ ! $2 ]; then
	IN_FILES=`ls *.sv.csv | tr '\n' ' '`
else
	IN_FILES=${*:2}
fi

echo "python /home/dennis.kao/tools/crg/crg.intersect_sv_reports.py -exon_bed=/home/dennis.kao/gene_panels/protein_coding_genes.exons.fixed.sorted.bed -o=${OUT} -i ${IN_FILES}"
python /home/dennis.kao/tools/crg/crg.intersect_sv_reports.py -exon_bed=/home/dennis.kao/gene_panels/protein_coding_genes.exons.fixed.sorted.bed -o=${OUT} -i ${IN_FILES}
