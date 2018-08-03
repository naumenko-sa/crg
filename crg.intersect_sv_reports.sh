#!/bin/bash

ls *.sv.csv | sed s/.sv.csv// > samples.txt

for sample in `cat samples.txt`
do
    cat $sample.sv.csv | sed 1d | awk -F '","' '{if ($2 == $6 ) print $0 }' > $sample.ins
    #head -n1 $sample.sv.csv > $sample.non_ins
    cat $sample.sv.csv | sed 1d | awk -F '","' '{if ($2 != $6 ) print $0 }' > $sample.non_ins
    
    cat $sample.non_ins | sed 1d | awk -F '","' -v smpl=$sample '{print $1"\t"$2"\t"$6"\t"smpl}' | sed s/"\""// | sort -k1,1 -k2,2n > $sample.bed
    #head -n1 $sample.sv.csv | awk '{print "superindex,"$0}' > $sample.sv.csv.indexed
    cat $sample.non_ins | awk -F '","' '{print $1"-"$2"-"$6"\","$0}' |sed s/'"'// | sort > $sample.non_ins.indexed
done

#overlapping with the first sample
#the problem is that this overlap merges all intervals, i.e. eats insertions
bedtools multiinter -header -i `ls 180_*bed` | awk '{print $1"-"$2"-"$3"\t"$0"\""}' | sed s/"\t"/'\",\"'/g | sed 1d | sort > 180.sv.mapped_to_first.csv

join -t '","' 180.sv.mapped_to_first.csv 180_120941Q.non_ins.indexed > test.csv


