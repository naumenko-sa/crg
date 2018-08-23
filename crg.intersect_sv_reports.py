import os
import csv
import re
import sys
import argparse

from datetime import datetime
from pybedtools import BedTool

csv.field_size_limit(sys.maxsize)

class StructuralVariant:

	'''Simple class to hold column values'''

	def __init__(self):

		#Standard BED fields
		self.chr = ""
		self.start = ""
		self.end = ""
		self.gene = ""

		#Custom implementation fields
		self.exons_spanned = ""

		#SVSCORES fields
		self.svlen = ""
		self.svtype = ""
		self.svmax = ""
		self.svsum = ""
		self.svtop5 = ""
		self.svtop10 = ""
		self.svmean = ""

		#TCAG Frequency
		self.dgv = ""

		#DECIPHER Population values
		self.del_obs = ""
		self.del_freq = ""
		self.del_se = ""
		self.dup_obs = ""
		self.dup_freq = ""
		self.dup_sta = ""
		self.obs = ""
		self.freq = ""
		self.se = ""
		self.cnv_type = ""
		self.samp_size = ""
		self.study = ""

		#DECIPHER CNVS values

	def make_decipher_link(self):
		return '=hyperlink("https://decipher.sanger.ac.uk/browser#q/{}:{}-{}")'.format(self.chr, self.start, self.end)

	def key(self):
		return (self.chr, self.start, self.end)

	def interval(self):
		return '{}:{}-{}:{}'.format(self.chr, self.start, self.end, self.svtype)

class StructuralVariantRecords:

	'''
		Holds groupings of StructuralVariant.

		Groupings in the dict grouped_sv follows this structure:
			grouped_sv[ref_interval][samp_name]

		Additional information about the ref_interval used to group other sample intervals is stored in:
			all_ref_interval_data[(chr, start, end)]
	'''

	def __init__(self, _sample_list):
		self.grouped_sv = {}
		self.all_ref_interval_data = {}
		self.sample_list = _sample_list

	def add_interval(self, ref_interval, new_interval, column_data, samp_name):
		'''
			Adds new_interval to the dict under the key ref_interval
			Column data relating to a ref_interval is stored as well
		'''
		if ref_interval not in self.grouped_sv.iterkeys():
			self.grouped_sv[ref_interval] = {}
			self.all_ref_interval_data[ref_interval] = column_data[ref_interval]

		if samp_name not in self.grouped_sv[ref_interval]:
			self.grouped_sv[ref_interval][samp_name] = []

		self.grouped_sv[ref_interval][samp_name].append(column_data[new_interval])

	def make_sample_list_index(self, interval):
		'''
			Makes string for LIST column
			e.g. 1;2;3;4;5;6;7
		'''
		index = []
		isthere_array = []

		for i, val in enumerate(self.sample_list):
			if val in interval:
				index.append(str(i+1))
				isthere_array.append("1")
			else:
				isthere_array.append("0")

		return (";".join(index), isthere_array)

	def make_sample_details(self, interval):
		'''
			Makes list for SAMPLENAME_details column
			e.g. 1:10334731-10334817:DEL;1:10334769-10334833:DUP
		'''
		all_samp_details = []

		for s in self.sample_list:

			samp_details = []

			if s in interval:
				for variant in interval[s]:
					samp_details.append(variant.interval())
				all_samp_details.append(';'.join(samp_details))
			else:
				all_samp_details.append("NA")

		return all_samp_details

	def make_header(self):
		fields = ["CHR", "START", "END", "N_SAMPLES", "LIST", "GENES", "LONGEST_SVTYPE", "SVLEN", "SVSCORE_MAX", "SVSCORE_SUM", "SVSCORE_TOP5", "SVSCORE_TOP10", "SVSCORE_MEAN", "DGV", "EXONS_SPANNED", "DECIPHER_LINK"]
		header = ",".join(fields)

		for s in self.sample_list:
			header = header + "," + s

		for s in self.sample_list:
			header = header + "," + "{}_details".format(s)

		return header + "\n"

	def get_longest_svtype(self, samples_in_interval):
		'''
			Determine the longest structural variant, then return its annotation
		'''
		longest = -1	#will always be overwritten on first loop iteration since we are using abs()
		longest_svtype = ""

		for samp in samples_in_interval:
			for variant in samples_in_interval[samp]:

				length = abs(int(variant.svlen))

				if length > longest:
					longest = length
					longest_svtype = variant.svtype

		return longest_svtype

	def BedTool_make_all_sv(self, bed_name):
		'''
			Creates a bed file containing all reference intervals and returns a BedTool object
			This is necessary since bedtools only operates on files
		'''
		with open(bed_name, "w") as f:
			for sv in self.grouped_sv:
				f.write('{}\t{}\t{}\t{}.\n'.format(sv[0], sv[1], sv[2], bed_name))

		return BedTool(bed_name)

	def calc_exons_spanned(self, exon_bed):
		'''
			Populates the field: exons_spanned for all reference intervals
		'''
		tmp_bed_name = "tmp_interval.bed"
		tmp_all_sv_bed_name = "tmp_all_sv.bed"

		exon_ref = BedTool(exon_bed)
		sample = self.BedTool_make_all_sv(tmp_all_sv_bed_name)

		for interval in sample:

			chr, start, end, gene = interval

			# create a temp bed file with 1 line - the interval of interest
			with open(tmp_bed_name, "w") as f:
				f.write('{}\t{}\t{}\n'.format(chr, start, end))

			tmp_bed = BedTool(tmp_bed_name)
			self.all_ref_interval_data[(chr, start, end)].exons_spanned = str(tmp_bed.intersect(exon_ref).count())

		os.remove(tmp_bed_name)
		os.remove(tmp_all_sv_bed_name)

	def decipher_pop_cnv_freq(decipher_population):
		
		cnv_freq = {}

		tmp_all_sv_bed_name = "tmp_all_sv.bed"
		decipher_pop_bed_name = "decipher_pop_cnv.bed"

		with open(decipher_population) as f:
			for sv in f:
				cnv_id, chr, start, end, del_obs, del_freq, del_se, dup_obs, dup_freq, dup_sta, obs, freq, se, cnv_type, samp_size, study = sv.split("\t")
				cnv_freq[(chr, start, end)] = [cnv_id, del_obs, del_freq, del_se, dup_obs, dup_freq, dup_sta, obs, freq, se, cnv_type, samp_size, study]

		with open(decipher_pop_bed_name, "w") as out_f:
			for interval in cnv_freq.iterkeys():
				out_f.write("\t".join(interval))

		all_sv_bed = self.BedTool_make_all_sv(tmp_all_sv_bed_name)
		pop_bed = BedTool(decipher_pop_bed_name)

	def write_results(self, outfile_name):
		'''
			A lot of string manipulation to generate the final CSV file line by line
		'''
		with open(outfile_name, "w") as out:

			out.write(self.make_header())

			for key in sorted(self.grouped_sv.iterkeys()):

				chr, start, end = key
				ref_interval = self.all_ref_interval_data[key]

				index, isthere = self.make_sample_list_index(self.grouped_sv[key])
				nsamples = str(len(self.grouped_sv[key]))
				svtype = self.get_longest_svtype(self.grouped_sv[key])
				samp_details = self.make_sample_details(self.grouped_sv[key])

				out_line = '{},{},{}\n'.format(",".join([str(chr), str(start), str(end), nsamples, index, ref_interval.gene, svtype, ref_interval.svlen, ref_interval.svmax, ref_interval.svsum, ref_interval.svtop5, ref_interval.svtop10, ref_interval.svmean, ref_interval.dgv, ref_interval.exons_spanned, ref_interval.make_decipher_link()]), ",".join(isthere), ",".join(samp_details))
				out.write(out_line)

def group_sv(a_bed_name, b_bed_name, column_data, all_sv_records):
	
	'''
		Groups intervals which overlap with a "reference" interval in sample a.

		Overlapping is defined as: each interval overlaps each other by >= 50% of their respective length.
		This gaurentees that grouping occurs when intervals are relatively the same size and overlap one another by a "significant" amount.

		Returns a BedTool instance containing all the intervals which did not meet the overlapping criteria.
	'''

	a_bed = BedTool(a_bed_name)
	b_bed = BedTool(b_bed_name)
	already_grouped_intervals = []

	# Find all SVs in b_bed which overlap a SV in a_bed by >=50% and are overlapped by >=50% by the same SV in first_sample_bed
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
		Taken from Sergey's previous script - crg.sv.merge_family.sh
		Converts each sample csv to a bed file using bash commands

		Returns list of sample file names with extension removed
	'''
	sample_list = []

	for s in input_files:

		if s.endswith('.sv.csv'):
			s=s[:-7]

		sample_list.append(s)
		cmd = "cat {}.sv.csv | sed 1d | awk -F \'\",\"\' -v smpl={} \'{{print $1\"\\t\"$2\"\\t\"$6\"\\t\"smpl}}\' | sed s/\"\\\"\"// | sort -k1,1 -k2,2n > {}.bed".format(s, s, s)
		os.popen(cmd)

	return sample_list

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

			chr, start, gt, svtype, svlen, end, sources, nsvt, gene, ann, svmax, svsum, svtop5, svtop10, svmean, dgv = line			
			key = (chr, start, end)

			if key not in column_data:
				newSV = StructuralVariant()

				newSV.chr = chr
				newSV.start = start
				newSV.svtype = svtype
				newSV.svlen = svlen
				newSV.end = end
				newSV.gene = gene
				newSV.svsum = svsum
				newSV.svmax = svmax
				newSV.svtop5 = svtop5
				newSV.svtop10 = svtop10
				newSV.svmean = svmean
				newSV.dgv = dgv

				column_data[key] = newSV

def main(exon_bed, outfile_name, input_files, decipher_population, decipher_pathogenic):
	'''
		Loop over input files performing bedtools intersect to get a grouping of intervals.
		Intervals which do not meet overlapping criteria are stored in a bed file in a pass_# folder for
		the next loop iteration. In my experience, these files only contain < 10 lines/intervals.

		See group_sv() for more details on what is defined as "overlapping criteria" 
	'''
	sample_list = csv2bed(input_files)
	all_sv_records = StructuralVariantRecords(sample_list)
	all_column_data = {}

	for i, s in enumerate(sample_list):

		current_dir = ""
		next_dir = ""

		if i != 0:
			current_dir = 'pass_{}/'.format(str(i))

		if i != len(sample_list)-1: # if not on last sample, create folder for input files in next processing pass
			next_dir = "pass_{}/".format(str(i+1))

			try:
				os.mkdir(next_dir)
			except OSError:
				pass # if directory already exists, just use it

		a_bed = current_dir + sample_list[i] + ".bed"
		a_csv = sample_list[i] + ".sv.csv"

		for j in range(i, len(sample_list)):

			b_bed = current_dir + sample_list[j] + ".bed"
			b_csv = sample_list[j] + ".sv.csv"
			parse_csv(b_csv, all_column_data)

			# print("pass {}: {} {}".format(str(i), a_bed, b_bed) )

			leftover_sv = group_sv(a_bed, b_bed, all_column_data, all_sv_records)

			if i != j:
				new_bed = next_dir + sample_list[j] + ".bed"
				make_bed_file(leftover_sv, new_bed) # store all leftover_sv in tmp dir for processing in next pass

	all_sv_records.calc_exons_spanned(exon_bed)
	all_sv_records.write_results(outfile_name)

if __name__ == '__main__':
	'''
		crg.intersect.sv.reports.py
		August 2018
		Center for Computational Medicine, SickKids, Dennis Kao

		Purpose: 
				To consolidate interval information across several samples and allow for analysis across a family/small population. In my work, intervals belong to the
				regions of structural variants. 
		How to use: 
				Navigate to a folder containing multiple .sv.csv files and launch the script. Then use Excel or some other spreadsheet software to filter the results.
		Output: 
				CSV file containing grouped intervals and metadata.
		Limitations: 
				Script is heavily IO bound because it relies on bedtools, a program which only operates on files. calc_exons_spanned() creates a file for each interval
				being calculated and thus is the most costly function. The main() loop creates files in pass_# folders containing intervals that did not meet the overlapping criteria
				but this doesn't really hinder performance since almost all intervals do tend to meet the overlapping criteria.

				Actual run time is acceptable for my uses. Researchers looking to run this code on a large sample size (>50) may want to look in to using RAMdisk.
		Typical runtime:
				<30 s for 2 samples
				>3 mins for 4 - 8 samples
		
		python crg.intersect_sv_reports.py -exon_bed="protein_coding_genes.exons.bed" -o="180.sv.family.csv" -i 180_123.sv.csv 180_444.sv.csv
	'''
	parser = argparse.ArgumentParser(description='Generates a structural variant report in CSV format for clincal research')
	parser.add_argument('-exon_bed', default="/home/naumenko/Desktop/reference_tables/protein_coding_genes.exons.fixed.sorted.bed", help='BED file containing fixed exon positions', required=True)
	parser.add_argument('-decipher_population', help='TSV file containing population CNVs frequencies from Decipher\'s internal database')
	parser.add_argument('-decipher_pathogenic', help='TSV file containing deleterious CNVs data from Decipher\'s internal database')
	parser.add_argument('-o', help='Output file name e.g. -o 180.sv.family.csv', required=True, type=str)
	parser.add_argument('-i', nargs='+', help='Input file names including .sv.csv extension, e.g. -i 180_230.sv.csv 180_231.sv.csv', required=True)
	args = parser.parse_args()

	print('crg.sv.merge_family.py started processing on ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))
	main(args.exon_bed, args.o, args.i, args.decipher_population, args.decipher_pathogenic)
	print('crg.sv.merge_family.py finished on ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))
