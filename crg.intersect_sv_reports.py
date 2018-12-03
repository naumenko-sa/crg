import os
import csv
import sys
import argparse
import shutil
import re
from datetime import datetime
from pybedtools import BedTool
from structuralvariant import GeneAnnotations, StructuralVariant, StructuralVariantRecords
csv.field_size_limit(sys.maxsize)

'''
	crg.intersect.sv.reports.py
	August 2018
	Center for Computational Medicine, SickKids, Dennis Kao

	Purpose: 
			To group and annotate SV's (structural variants) of similar location and size across several samples.
			Enables easier spreadsheet analysis across a cohort.
	How to use:
			Navigate to a folder containing multiple .sv.csv file and launch the script. Then use Excel or some other spreadsheet software to filter the results.
	Output:
			CSV file containing grouped intervals and annotated columns.
	Typical runtime on 2012 Macbook Pro i5 3210M, 8GB DDR3-1600:
				<30s
	Requirements:
			Python 3
			pybedtools
			sqlite3
	To run:
			python crg.intersect_sv_reports.py -exon_bed="protein_coding_genes.exons.bed" -o="180.sv.family.csv" -i 180_123.sv.csv 180_444.sv.csv
'''

def group_sv(a_bed, b_bed, column_data, all_sv_records):
	'''
		Groups SV's which overlap with a "reference" interval in sample a.

		Overlapping is defined as: a recipricol overlap of >=50%.
		This gaurentees grouping occurs between SV's of similar size and position.

		Returns a BedTool instance containing all the intervals which did not meet the overlapping criteria.
	'''
	already_grouped_intervals = [] #handles edge case where intervals within the same sample overlap one another. This can lead to an SV appearing in more than 1 grouping.
	overlapping_sv = a_bed.intersect(b_bed, wa=True, wb=True, F=0.5, f=0.5)

	for l in overlapping_sv:
		chr, start, end, ref_name, samp_chr, samp_start, samp_end, samp_name = l

		ref_interval = (chr, start, end)
		new_interval = (samp_chr, samp_start, samp_end)

		if new_interval not in already_grouped_intervals:
			all_sv_records.add_interval(ref_interval, new_interval, column_data, samp_name)
			already_grouped_intervals.append(new_interval)

	return b_bed.intersect(a_bed, F=0.5, f=0.5, v=True)

def parse_sample_csv(sample_csv, column_data):
	'''
		Parse a sample's .sv.csv file and stores meta data associated with an interval to all_ref_interval_data.
		Returns a bedtools instance containing all intervals within a sample.
	'''
	intervals = []
	sample_name = sample_csv[:-7] if sample_csv.endswith(".sv.csv") else sample_csv
	with open(sample_csv) as f:
		try: # skip the header
			next(f)
		except ValueError: #file is empty
			return BedTool([])

		for line in csv.reader(f, delimiter=",", quotechar="\""):
			if line:

				chr, start, genotype, svtype, svlen, end, sources, nsvt, genes, ann, svmax, svsum, svtop5, svtop10, svmean, dgv = line	
				dgv = "0" if dgv == "NA" else dgv		
				key = (chr, start, end)
				intervals.append((chr, start, end, sample_name))
				
				if key not in column_data:
					newSV = StructuralVariant(chr, start, end, svtype, genotype, svlen, svsum, svmax, svtop5, svtop10, svmean, dgv)
					for gene_name in set(re.split('[&;]+', genes.strip())): #set ensures uniqueness - the same gene doesn't get added twice
						if gene_name:
							newSV.add_gene(gene_name)
					column_data[key] = newSV	

	return BedTool(intervals)

def main(exon_bed, hgmd_db, hpo, exac, omim, outfile_name, input_files):
	'''
		Loop over input files performing bedtools intersect to get a grouping of intervals.
		See group_sv() for more details on what is defined as "overlapping criteria".

		Then, annotate and create CSV report.
	'''
	all_sv_records = StructuralVariantRecords([ sample[:-7] if sample.endswith('.sv.csv') else sample for sample in input_files ])
	all_column_data = {}
	passes = [[parse_sample_csv(f, all_column_data) for f in input_files], ]

	print("Grouping like structural variants ...")
	for current_pass in passes:
		a_bedtool = current_pass[0]
		next_pass = []
		for b_bedtool in current_pass:
			leftover_sv = group_sv(a_bedtool, b_bedtool, all_column_data, all_sv_records)
			if leftover_sv: next_pass.append(leftover_sv)
		if next_pass: passes.append(next_pass) #keep on processing if there are ungrouped SVs

	print('Annotating structural variants ...')
	all_sv_records.annotate(exon_bed, hgmd_db, hpo, exac, omim)
	print('Writing results to file ...')
	all_sv_records.write_results(outfile_name)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Generates a structural variant report in CSV format for clincal research')
	parser.add_argument('-exon_bed', default="/home/dennis.kao/gene_panels/protein_coding_genes.exons.fixed.sorted.bed", help='BED file containing fixed exon positions', required=True)
	parser.add_argument('-hgmd', help='HGMD Pro database file', required=True, type=str)
	parser.add_argument('-hpo', help='Tab delimited file containing gene names and HPO terms', type=str)
	parser.add_argument('-exac', help='ExAC tab delimited file containing gene names and scores', type=str, required=True)
	parser.add_argument('-omim', help='OMIM tab delimited file containing gene names and scores', type=str, required=True)
	parser.add_argument('-o', help='Output file name e.g. -o 180.sv.family.csv', required=True, type=str)
	parser.add_argument('-i', nargs='+', help='Input file names including .sv.csv extension, e.g. -i 180_230.sv.csv 180_231.sv.csv', required=True)
	args = parser.parse_args()

	print('crg.intersect_sv_reports.py started processing on ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))
	main(args.exon_bed, args.hgmd, args.hpo, args.exac, args.omim, args.o, args.i)
	print('crg.intersect_sv_reports.py finished on ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))