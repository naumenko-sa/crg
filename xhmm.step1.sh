#!/bin/bash
#PBS -l walltime=23:59:59,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=50g,mem=50g

#http://atgu.mgh.harvard.edu/xhmm/tutorial.shtml#params_file

#bam.list
#EXOME.interval_list is a list of exon (could be limited only to 1 chr)
#chr:start-end

# GATK4 does not have DepthOfCoverage
# -U ALLOW_N_CIGAR_READS for RNA-seq

reference=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa
xhmm_home=/hpf/largeprojects/ccmbio/naumenko/tools/xhmm/statgen-xhmm-cc14e528d909

gatk3 -Xmx3072m -T DepthOfCoverage \
-I bam.$n.list -L $interval_list \
-R $reference \
-dt BY_SAMPLE -dcov 5000 -l INFO --omitDepthOutputAtEachBase --omitLocusTable \
--minBaseQuality 0 --minMappingQuality 20 --start 1 --stop 5000 --nBins 200 \
--includeRefNSites \
--countType COUNT_FRAGMENTS \
-o group.$n.DATA \
-U ALLOW_N_CIGAR_READS
