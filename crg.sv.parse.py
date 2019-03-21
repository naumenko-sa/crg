#!/bin/env python

import csv
import sys
import re

vcf = sys.argv[1]

dgv_dict = {}
if(len(sys.argv) > 2):
    dgv_report = sys.argv[2]
    with open(dgv_report) as f_dgv_report:
	for line in f_dgv_report:
	    if not line.startswith("SAMPLE_ID"):
		buf = line.strip()
		fields = buf.split('\t')
		key = fields[1]+'-'+fields[2]+'-'+fields[3]
		dgv_dict[key] = fields[47]
		#print(key+'\t'+dgv_dict[key])

colnames = ['CHR','POS','GT','SVTYPE','SVLEN','END','SOURCES','NUM_SVTOOLS','GENES','ANN','SVSCOREMAX','SVSCORESUM','SVSCORETOP5','SVSCORETOP10','SVSCOREMEAN','DGV']

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

	    genes = []
	    for info_value in info_fields:
		pair = info_value.split('=')
		if (len(pair) == 2):
			info_dict[pair[0]] = pair[1]
			if (pair[0] == 'ANN'):
			    ann_buf = pair[1].split(',')
			    for ann in ann_buf:
				ann_buf1 = ann.split('|')
				genes.append(ann_buf1[3])
	    info_dict['GENES'] = ';'.join(set(genes))
				
	    for field in colnames[3:15]:
		if field in info_dict:
		    values.append(info_dict[field])
		else:
		    values.append('NA')

	    key = values[0]+'-'+values[1]+'-'+info_dict['END']
	    #print(key+"\t"+dgv_dict[key])
	    if dgv_dict:
		if key in dgv_dict:
		    info_dict['DGV'] = dgv_dict[key]
    		else:
		    info_dict['DGV'] = '0'
	    else:
		info_dict['DGV'] = '0'

	    if len(sys.argv) > 2:
		values.append(info_dict['DGV'])
	    else:
		values.append('NA')

	    print '"'+'","'.join(values)+'"'

f_vcf.close()

