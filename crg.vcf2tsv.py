#!/bin/env python2

#converts vcf to tsv to send to TCAG for frequency annotation
# TCAG wants:
# Sample ID
# Chromosome
# Start
# End
# SV type  (they work only with deletion and insertions)
# (additional columns..)

import csv
import sys
import re
import vcf

#vcf = sys.argv[1]

vcf_reader = vcf.Reader(open(sys.argv[1]),'r')
sample_id = vcf_reader.samples[0]

colnames = ['SAMPLE_ID','CHR','POS','END','SVTYPE']

print '\t'.join(colnames)

with open(sys.argv[1]) as f_vcf:
    for line in f_vcf:
	if not line.startswith('#'):
	    buf = line.strip()
	    vcf_fields = buf.split()
	    values = list()
	    values.append(sample_id)
	    values.append(vcf_fields[0])
	    values.append(vcf_fields[1])
	    
	    info_string = vcf_fields[7]
	    
	    info_fields = info_string.split(';')
	    info_dict = {}
	    for info_value in info_fields:
		pair = info_value.split('=')
		if(len(pair)==2):
		    info_dict[pair[0]] = pair[1]
	
	    for field in colnames[3:]:
		values.append(info_dict[field])

	    print '\t'.join(values)

f_vcf.close()
