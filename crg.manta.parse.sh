#!/bin/bash

#filter SVs produced by Manta (default Illumina pipeline)

#for f in *.vcf.gz;do bcftools view -f "PASS,." -R sle.bed $f | bcftools sort - | bgzip > `echo $f | sed s/vcf.gz/sorted.vcf.gz/`;done;

#for f in `cat samples.txt`;do vcfanno gene_name.toml ${f}_S1.SV.sorted.vcf.gz | bgzip > $f.sv.vcf.gz; tabix $f.sv.vcf.gz;done;

for f in `cat samples.txt`
do 
    bcftools query -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t%REF\t%ALT\t%gene\t[%GT]\n' $f.sv.vcf.gz | awk -v smpl=$f '{print smpl"\t"$0}' > $f.tsv
done

rm SV.tsv
echo -e 'SAMPLE\tCHROM\tPOS\tEND\tSVLEN\tSVTYPE\tREF\tALT\tgene\tGT' > SV.ts
cat *.tsv >> SV.ts
mv SV.ts SV.tsv
