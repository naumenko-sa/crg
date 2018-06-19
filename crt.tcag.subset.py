#!/bin/env python

# subset tcag report using a bed file

import csv
import sys
import re

original_report = sys.argv[1]

gene_panel_file_name = sys.argv[2]
bed = []

with open(gene_panel_file_name,'rb') as f_gene_panel:
    for line in f_gene_panel:
	record = line.strip().split('\t')
	bed.append(record)
f_gene_panel.close()
	
with open(original_report,'rb') as f_original_report:
    for line in f_original_report:
	if line.startswith('#'):
	    header_line = line.strip().replace('\t',',')
	    print(header_line)
	else:
	    record = line.strip().split('\t')
	    chrom = record[0].replace('chr','')
	    start = record[1]
	    for interval in bed:
		if ((interval[0] == chrom) and (start>=interval[1]) and (start<=interval[2])):
		    print '"'+'","'.join(record)+'"'

f_original_report.close()
