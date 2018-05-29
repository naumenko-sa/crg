#!/bin/env python2

import csv
import sys
import re

vcf = sys.argv[1]

colnames = ['CHR','POS','GT','SVTYPE','SVLEN','END','SOURCES','NUM_SVTOOLS','ANN']
#add genes - parse ANN

print ','.join(colnames)

with open(vcf) as f_vcf:
    for line in f_vcf:
	if not line.startswith('#'):
	    buf = line.strip()
	    vcf_fields = buf.split()
	    values = list()
	    values.append(vcf_fields[0])
	    values.append(vcf_fields[1])
	    
	    if (vcf_fields[9] == '0/1'):
		genotype = 'HET'
	    else:
		genotype = 'HOM'
	    values.append(genotype)
	    
	    info_string = vcf_fields[7]
	    
	    info_fields = info_string.split(';')
	    info_dict = {}
	    info_dict['ANN']='NA'
	    for info_value in info_fields:
		pair = info_value.split('=')
		if (len(pair) == 2):
			info_dict[pair[0]] = pair[1]
	
	    for field in colnames[3:]:
		values.append(info_dict[field])

	    print '"'+'","'.join(values)+'"'

f_vcf.close()
