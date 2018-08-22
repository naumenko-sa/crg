import os
import csv
import re
import sys
import argparse

from datetime import datetime
from pybedtools import BedTool

csv.field_size_limit(sys.maxsize)

class StructuralVariantRecords:

	def __init__(self, _sample_list):
		self.all_sv = {}
		self.all_sv_column_data = {}
		self.sample_list = _sample_list

	def add_interval(self, ref_interval, new_interval, new_svlen, new_svtype, samp_name):
		if ref_interval not in self.all_sv.iterkeys():
			self.all_sv[ref_interval] = {}

		if samp_name not in self.all_sv[ref_interval]:
			self.all_sv[ref_interval][samp_name] = []

		self.all_sv[ref_interval][samp_name].append((new_interval, new_svlen, new_svtype))

	def parse_column_data(self, sample_csv):
		'''
			Parse a sample's CSV and store column data associated with an interval to all_sv_column_data.
		'''
		with open(sample_csv) as f:

			try: # skip the header
				next(f)
			except ValueError: #file is empty
				return

			for line in csv.reader(f, delimiter=",", quotechar="\""):

				if not line:
					continue

				chr, start, gt, svtype, svlen, end, sources, nsvt, gene, ann, svmax, svsum, svtop5, svtop10, svmean, dvg = line
				key = (chr, start, end)

				if key not in self.all_sv_column_data:
					self.all_sv_column_data[key] = (svlen, svmax, svsum, svtop5, svtop10, svmean, dvg, gene)

	def make_sample_list_index(self, interval):
		'''
			LIST column
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
			SAMPLENAME_details column
			e.g. 1:10334731-10334817:DEL;1:10334769-10334833:DUP
		'''

		#print(interval)

		all_samp_details = []

		for s in self.sample_list:

			samp_details = []

			if s in interval:
				for variant in interval[s]:
					chr, start, end = variant[0]
					svtype = variant[2]

					samp_details.append('{}:{}-{}:{}'.format(chr, start, end, svtype))
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
				chr, start, end = variant[0]
				svlen, svtype = variant[1], variant[2]

				length = abs(int(svlen))

				if length > longest:
					longest = length
					longest_svtype = svtype

		return longest_svtype

	def make_spanning_exon_dict(self, exon_bed):
		'''
			Returns a dict where each interval maps to the number of exons it spans
		'''
		tmp_bed_name = "tmp_interval.bed"
		tmp_all_sv_bed_name = "tmp_all_sv.bed"
		overlapping_exon_count = {}

		# generate temporary bed file containing all intervals
		with open(tmp_all_sv_bed_name, "w") as f:
			for sv in self.all_sv:
				f.write('{}\t{}\t{}\t.\n'.format(sv[0], sv[1], sv[2]))

		exon_ref = BedTool(exon_bed)
		sample = BedTool(tmp_all_sv_bed_name)

		for interval in sample:

			chr, start, end, gene = interval

			# create a temp bed file with 1 line - the interval of interest
			with open(tmp_bed_name, "w") as f:
				f.write('{}\t{}\t{}\n'.format(interval.chrom, interval.start, interval.end))

			tmp_bed = BedTool(tmp_bed_name)
			overlapping_exon_count[(chr, start, end)] = tmp_bed.intersect(exon_ref).count()

		os.remove(tmp_bed_name)
		os.remove(tmp_all_sv_bed_name)

		return overlapping_exon_count

	def write_results(self, outfile_name, exon_bed):
		'''
			A lot of string manipulation to generate the final CSV file line by line
		'''
		overlapping_exon_count = self.make_spanning_exon_dict(exon_bed)

		with open(outfile_name, "w") as out:

			out.write(self.make_header())

			for key in sorted(self.all_sv.iterkeys()):

				chr, start, end = key

				index, isthere = self.make_sample_list_index(self.all_sv[key])
				nsamples = str(len(self.all_sv[key]))
				svtype = self.get_longest_svtype(self.all_sv[key])
				samp_details = self.make_sample_details(self.all_sv[key])
				decipher_link = '=hyperlink("https://decipher.sanger.ac.uk/browser#q/{}:{}-{}")'.format(chr, start, end)
				n_exon_spanned = str(overlapping_exon_count[key])
				svlen, svmax, svsum, svtop5, svtop10, svmean, dvg, gene = self.all_sv_column_data[key]

				out_line = '{},{},{}\n'.format(",".join([str(chr), str(start), str(end), nsamples, index, gene, svtype, svlen, svmax, svsum, svtop5, svtop10, svmean, dvg, n_exon_spanned, decipher_link]), ",".join(isthere), ",".join(samp_details))
				out.write(out_line)

def determine_intersect(first_sample_bed, other_sample_bed, other_sample_sv_info, all_sv_records):
	
	'''
		Populate all_sv dict with overlapping SVs from another bed file.

		Returns intervals which did not meet overlapping criteria.
	'''

	a_bed = BedTool(first_sample_bed)
	b_bed = BedTool(other_sample_bed)
	grouped_intervals = []

	# find all SVs in b_bed which overlap a SV in a_bed by >=50% and are overlapped by >=50% by the same SV in first_sample_bed
	# bedtools doc for more param info: https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html?highlight=intersect
	overlapping_sv = a_bed.intersect(b_bed, wa=True, wb=True, F=0.5, f=0.5)

	for l in overlapping_sv:

		chr, start, end, ref_name, samp_chr, samp_start, samp_end, samp_name = l

		ref_interval = (chr, start, end)
		new_interval = (samp_chr, samp_start, samp_end)

		if new_interval not in grouped_intervals:
			all_sv_records.add_interval(ref_interval, new_interval, other_sample_sv_info[new_interval][0], other_sample_sv_info[new_interval][1], samp_name)
			grouped_intervals.append(new_interval)

	# report all SV in B which did not meet intersection criteria above
	return b_bed.intersect(a_bed, F=0.5, f=0.5, v=True)

def make_svtype_svlen_dict(scsv):

	'''
		Store each interval's SVTYPE from CSV file.
		This info is later used in get_longest_svtype().

		vdict format: vdict[(chr, start, end)] = (svlen, svtype)
	'''

	vdict = {}

	with open(scsv) as f:

		try: # skip the header
			next(f)
		except ValueError: # file is empty
			return vdict

		for line in csv.reader(f, delimiter=",", quotechar="\""):

			if not line:
				continue

			chr, start, gt, svtype, svlen, end, sources, nsvt, genes, ann, svmax, svsum, svtop5, svtop10, svmean, dvg = line

			key = (chr, start, end)
			vdict[key] = (svlen, svtype)

	return vdict

def csv2bed(input_files):
	'''
		Taken from sergey's previous script - crg.sv.merge_family.sh
		Converts each sample csv to a bed file
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

def make_bed_file(interval, bed_dir):
	with open(bed_dir, "w") as f:
		for i in interval:
			f.write(str(i))

def main(exon_bed, outfile_name, input_files, decipher_population, decipher_pathogenic):

	sample_list = csv2bed(input_files)
	all_sv_records = StructuralVariantRecords(sample_list)

	for i, s in enumerate(sample_list):

		current_dir = ""
		next_dir = ""

		if i != 0:
			current_dir = 'pass_{}/'.format(str(i))

		if i != len(sample_list)-1: # if not on last sample, create folder for input files in next processing loop
			next_dir = "pass_{}/".format(str(i+1))

			try:
				os.mkdir(next_dir)
			except OSError:
				pass # if directory already exists, just use it

		a_bed = current_dir + sample_list[i] + ".bed"
		a_csv = sample_list[i] + ".sv.csv"

		all_sv_records.parse_column_data(a_csv)

		for j in range(i, len(sample_list)):

			b_bed = current_dir + sample_list[j] + ".bed"
			b_csv = sample_list[j] + ".sv.csv"
			b_svtype_svlen_dict = make_svtype_svlen_dict(b_csv)

			# print("pass {}: {} {}".format(str(i), a_bed, b_bed) )

			leftover_sv = determine_intersect(a_bed, b_bed, b_svtype_svlen_dict, all_sv_records)

			if i != j:
				new_bed = next_dir + sample_list[j] + ".bed"
				make_bed_file(leftover_sv, new_bed) # store all leftover_sv in tmp dir for processing in next pass

	all_sv_records.write_results(outfile_name, exon_bed)

if __name__ == '__main__':
	'''
		How to use: Navigate to a folder containing multiple sv.csv files and launch the script.
		Program generates a CSV file named after outfile variable.

		python crg.intersect_sv_reports.py -exon_bed="protein_coding_genes.exons.bed" -o="180.sv.family.csv" -i 180_123.sv.csv 180_444.sv.csv
	'''
	parser = argparse.ArgumentParser(description='Generates a structural variant report in CSV format for clincal research')
	parser.add_argument('-exon_bed', default="/home/naumenko/Desktop/reference_tables/protein_coding_genes.exons.fixed.sorted.bed", help='BED file containing fixed exon positions', required=True)
	parser.add_argument('-decipher_population', help='TSV file containing population CNVs frequencies from Decipher\'s internal database')
	parser.add_argument('-decipher_pathogenic', help='TSV file containing deleterious CNVs descriptions from Decipher\'s internal database')
	parser.add_argument('-o', help='Output file name e.g. -o 180.sv.family.csv', required=True, type=str)
	parser.add_argument('-i', nargs='+', help='Input file names including .sv.csv extension, e.g. -i 180_230.sv.csv 180_231.sv.csv', required=True)
	args = parser.parse_args()

	print('crg.sv.merge_family.py started processing on ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))
	main(args.exon_bed, args.o, args.i, args.decipher_population, args.decipher_pathogenic)
	print('crg.sv.merge_family.py finished on ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))
