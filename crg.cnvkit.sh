#!/bin/bash

#cnvkit variant calling based on bcbio log

#PBS -l walltime=10:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=21g,mem=21g

#vpath=/hpf/largeprojects/ccmbio/naumenko/validation/2018-03-01_venter_SV
#sample=912R_A337376

#cnvkit.py access /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa -o access.grch37.bed

#export TMPDIR=. && cnvkit.py target $target --split \
#    -o $sample.target.bed --avg-size 250
#export TMPDIR=. && cnvkit.py antitarget $target -g access.grch37.bed \
#    -o $sample.antitarget.bed --avg-size 1000000
export TMPDIR=. && cnvkit.py reference -f /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa \
    -o FlatReference.cnn \
    -t $sample.target.bed \
    -a $sample.antitarget.bed

#export TMPDIR=. && cnvkit.py fix -o 912R_A337376-normalized.cnr $vpath/work/structural/912R_A337376/bins/912R_A337376-target-coverage.cnn \
#    $vpath/work/structural/912R_A337376/bins/912R_A337376-antitarget-coverage.cnn \
#    background-0-cnvkit.cnn --sample-id 912R_A337376
#unset R_HOME && unset R_LIBS && export PATH=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/anaconda/bin:$PATH && export TMPDIR=. && cnvkit.py segment -p 7 \
#    -o 912R_A337376.cns 912R_A337376-normalized.cnr --threshold 0.00001
#cnvkit.py export seg -o 912R_A337376.seg 912R_A337376.cns
#cnvkit.py gainloss -s 912R_A337376.cns -o 912R_A337376-gainloss.txt 912R_A337376-normalized.cnr
#cnvkit.py segmetrics --median --iqr --ci --pi -s 912R_A337376.cns -o 912R_A337376-segmetrics.txt 912R_A337376-normalized.cnr --alpha 0.1 --bootstrap 50
#cnvkit.py call --filter cn --ploidy 2 -o 912R_A337376-call.cns 912R_A337376.cns
#cnvkit.py export bed --sample-id 912R_A337376 --ploidy 2 -o 912R_A337376-call.bed 912R_A337376-call.cns
#cnvkit.py export vcf --sample-id 912R_A337376 --ploidy 2 -o 912R_A337376-call.vcf 912R_A337376-call.cns
#cat 912R_A337376-call.bed | grep -v ^track | grep -v ^browser | grep -v "^#" | sort -V -T . -k1,1 -k2,2n | \
#    bedtools closest -g <(cut -f1,2 /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa.fai | \
#    sort -V -T . -k1,1 -k2,2n) -d -t all -a - -b <(sort -V -T . -k1,1 -k2,2n /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/rnaseq/ref-transcripts.bed) -t first | awk -F$'\t' -v OFS='\t' '{if ($NF > 10000) $9 = "."} {print}' | cut -f 1-10 | \
#    bedtools merge -i - -c 4,5,9 -o distinct,distinct,distinct -delim ',' -d -10 > 912R_A337376-call-annotated.bed
#unset JAVA_HOME && export PATH=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/anaconda/bin:$PATH &&  /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/anaconda/bin/snpEff -Xms750m -Xmx20g \
#-Djava.io.tmpdir=. eff -dataDir /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/snpeff -hgvs -noLog -i vcf -o vcf -csvStats \
#912R_A337376-call-effects-stats.csv -s 912R_A337376-call-effects-stats.html GRCh37.75 912R_A337376-call.vcf  > 912R_A337376-call-effects.vcf
#export TMPDIR=/hpf/largeprojects/ccmbio/naumenko/validation/venter/work/bcbiotx/tmpoCV55r && /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/anaconda/bin/cnvkit.py diagram -s /hpf/largeprojects/ccmbio/naumenko/validation/venter/work/structural/912R_A337376/cnvkit/912R_A337376-chromfilter.cns -o /hpf/largeprojects/ccmbio/naumenko/validation/venter/work/bcbiotx/tmpoCV55r/912R_A337376-normalized-diagram.pdf /hpf/largeprojects/ccmbio/naumenko/validation/venter/work/structural/912R_A337376/bins/912R_A337376-normalized-chromfilter.cnr
#export TMPDIR=/hpf/largeprojects/ccmbio/naumenko/validation/venter/work/bcbiotx/tmpeL4KlG && /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/anaconda/bin/cnvkit.py scatter -s /hpf/largeprojects/ccmbio/naumenko/validation/venter/work/structural/912R_A337376/cnvkit/912R_A337376-chromfilter.cns -o /hpf/largeprojects/ccmbio/naumenko/validation/venter/work/bcbiotx/tmpeL4KlG/912R_A337376-normalized-scatter_global.pdf /hpf/largeprojects/ccmbio/naumenko/validation/venter/work/structural/912R_A337376/bins/912R_A337376-normalized-chromfilter.cnr
