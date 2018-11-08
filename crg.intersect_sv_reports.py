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

def group_sv(a_bed_name, b_bed_name, column_data, all_sv_records):
	
	'''
		Groups SV's which overlap with a "reference" interval in sample a.

		Overlapping is defined as: a recipricol overlap of >=50%.
		This gaurentees grouping occurs between SV's of similar size and position.

		Returns a BedTool instance containing all the intervals which did not meet the overlapping criteria.
	'''

	a_bed = BedTool(a_bed_name)
	b_bed = BedTool(b_bed_name)
	already_grouped_intervals = [] #list handles edge case where intervals within the same sample overlap one another. This can lead to an SV appearing in more than 1 grouping.

	overlapping_sv = a_bed.intersect(b_bed, wa=True, wb=True, F=0.5, f=0.5)

	for l in overlapping_sv:

		chr, start, end, ref_name, samp_chr, samp_start, samp_end, samp_name = l

		ref_interval = (chr, start, end)
		new_interval = (samp_chr, samp_start, samp_end)

		if new_interval not in already_grouped_intervals:
			all_sv_records.add_interval(ref_interval, new_interval, column_data, samp_name)
			already_grouped_intervals.append(new_interval)

	return b_bed.intersect(a_bed, F=0.5, f=0.5, v=True)

def csv2bed(input_files):
	'''
		Taken from Sergey's script - crg.sv.merge_family.sh
		Converts each sample csv to a bed file
	'''	
	for s in input_files:
		cmd = "cat {}.sv.csv | sed 1d | awk -F \'\",\"\' -v smpl={} \'{{print $1\"\\t\"$2\"\\t\"$6\"\\t\"smpl}}\' | sed s/\"\\\"\"// | sort -k1,1 -k2,2n > {}.bed".format(s, s, s)
		os.popen(cmd)

def make_bed_file(intervals, bed_dir):
	'''
		Write all intervals from a BedTools instance to a .bed file
	'''
	with open(bed_dir, "w") as f:
		for i in intervals:
			f.write(str(i))

def parse_csv(sample_csv, column_data):
	'''
		Parse a sample's CSV and stores column data associated with an interval to all_ref_interval_data.
	'''
	with open(sample_csv) as f:

		try: # skip the header
			next(f)
		except ValueError: #file is empty
			return

		for line in csv.reader(f, delimiter=",", quotechar="\""):
			if not line:
				continue

			chr, start, genotype, svtype, svlen, end, sources, nsvt, genes, ann, svmax, svsum, svtop5, svtop10, svmean, dgv = line			
			key = (chr, start, end)

			if key not in column_data:
				newSV = StructuralVariant(chr, start, end, svtype, genotype, svlen, svsum, svmax, svtop5, svtop10, svmean, dgv)
				for gene_name in set(re.split('[&;]+', genes.strip())): #set ensures uniqueness - the same gene doesn't get added twice
					if gene_name:
						newSV.add_gene(gene_name)
				column_data[key] = newSV

def main(exon_bed, hgmd_db, outfile_name, input_files):
	'''
		Loop over input files performing bedtools intersect to get a grouping of intervals.
		Intervals which do not meet overlapping criteria are stored in a bed file in a pass_# folder for
		the next loop iteration. In my experience, these files contain < 10 lines/intervals.

		See group_sv() for more details on what is defined as "overlapping criteria" 
	'''
	def cleanup(tmp_bed, tmp_dir):
		for f in tmp_bed:
			os.remove(f)

		for d in tmp_dir:
			shutil.rmtree(d)

	sample_list = [ sample[:-7] if sample.endswith('.sv.csv') else sample for sample in input_files ]
	all_sv_records = StructuralVariantRecords(sample_list)
	all_column_data = {}
	tmp_dir = []

	csv2bed(sample_list)

	for i, s in enumerate(sample_list):

		current_dir = ""
		next_dir = ""

		if i != 0:
			current_dir = 'pass_{}/'.format(str(i))

		if i != len(sample_list)-1: # if not on last sample, create folder for input files in next processing pass
			next_dir = "pass_{}/".format(str(i+1))
			tmp_dir.append(next_dir)

			try:
				os.mkdir(next_dir)
			except OSError:
				pass # if directory already exists, just use it

		a_bed = "%s%s.bed" % (current_dir, s)

		for j in range(i, len(sample_list)):

			b_bed = "%s%s.bed" % (current_dir, sample_list[j])
			b_csv = "%s.sv.csv" % sample_list[j]
			parse_csv(b_csv, all_column_data)

			# print("pass {}: {} {}".format(str(i), a_bed, b_bed) )

			leftover_sv = group_sv(a_bed, b_bed, all_column_data, all_sv_records)

			if i != j:
				new_bed = "%s%s.bed" % (next_dir, sample_list[j])
				make_bed_file(leftover_sv, new_bed) # store all leftover_sv in tmp dir for processing in next pass

	all_sv_records.annotate(exon_bed, hgmd_db)
	cleanup(['%s.bed' % s for s in sample_list], tmp_dir)
	all_sv_records.write_results(outfile_name)

if __name__ == '__main__':
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
		Limitations:
				calc_exons_spanned() creates a file for each SV being calculated. Increases run-time by a minute or more.
		Typical runtime:
				<30 s for <=2 samples
				>3 mins for 4 - 8 samples
		Requirements:
			python 3.7 - needed since dictionary order is gaurenteed
			pybedtools
			sqlite3
			unix file system
		
		python crg.intersect_sv_reports.py -exon_bed="protein_coding_genes.exons.bed" -o="180.sv.family.csv" -i 180_123.sv.csv 180_444.sv.csv
	'''
	parser = argparse.ArgumentParser(description='Generates a structural variant report in CSV format for clincal research')
	parser.add_argument('-exon_bed', default="/home/dennis.kao/gene_panels/protein_coding_genes.exons.fixed.sorted.bed", help='BED file containing fixed exon positions', required=True)
	parser.add_argument('-hgmd', help='Path to hgmd_pro.db file', required=True, type=str)
	parser.add_argument('-o', help='Output file name e.g. -o 180.sv.family.csv', required=True, type=str)
	parser.add_argument('-i', nargs='+', help='Input file names including .sv.csv extension, e.g. -i 180_230.sv.csv 180_231.sv.csv', required=True)
	args = parser.parse_args()

	print('crg.intersect_sv_reports.py started processing on ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))
	main(args.exon_bed, args.hgmd, args.o, args.i)
	print('crg.intersect_sv_reports.py finished on ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))
