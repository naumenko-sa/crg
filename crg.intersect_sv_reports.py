import os
import csv
import re
import sys
import argparse

from datetime import datetime
from pybedtools import BedTool

csv.field_size_limit(sys.maxsize)

def determine_intersect(first_sample_bed, other_sample_bed, other_sample_sv_info, all_sv):
	
	'''
		Populate all_sv dict with overlapping SVs from another bed file.

		Returns intervals which did not meet overlapping criteria.
	'''

	a_bed = BedTool(first_sample_bed)
	b_bed = BedTool(other_sample_bed)

	# find all SVs in other_sample_bed which overlap a SV in first_sample_bed by >=50% and are overlapped by >=50% by the same SV in first_sample_bed
	# bedtools doc for more param info: https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html?highlight=intersect
	overlapping_sv = a_bed.intersect(b_bed, wa=True, wb=True, F=0.5, f=0.5)

	for l in overlapping_sv:

		chr, start, end, ref_name, samp_chr, samp_start, samp_end, samp_name = l

		interval = (chr, start, end)
		other_interval = (samp_chr, samp_start, samp_end)

		if interval not in all_sv:
			all_sv[interval] = {}

		if samp_name not in all_sv[interval]:
			all_sv[interval][samp_name] = []

		all_sv[interval][samp_name].append((other_interval, other_sample_sv_info[other_interval]))

	# report all SV in B which did not meet intersection criteria above
	return b_bed.intersect(a_bed, F=0.5, f=0.5, v=True)

def save_sv_column_data(scsv, all_sv_info):

	'''
		Parse the first sample's CSV and store column data associated with interval.

		variant_info holds SV_SCOREs and other column info from the CSV file
	'''
	with open(scsv) as f:

		next(f)	# skip the header

		for line in csv.reader(f, delimiter=",", quotechar="\""):

			if not line:
				continue

			chr, start, gt, svtype, svlen, end, sources, nsvt, gene, ann, svmax, svsum, svtop5, svtop10, svmean, dvg = line
			key = (chr, start, end)

			if key not in all_sv_info:
				all_sv_info[key] = (svlen, svmax, svsum, svtop5, svtop10, svmean, dvg, gene)
			# elif key in all_sv_info:
			# 	if all_sv_info[key][0] 

def make_svtype_svlen_dict(scsv):

	'''
		Store each interval's SVTYPE from CSV file.
		This info is later used in get_longest_svtype().

		vdict format: vdict[(chr, start, end)] = svtype
	'''

	with open(scsv) as f:

		vdict = {}

		next(f)	#skip the header

		for line in csv.reader(f, delimiter=",", quotechar="\""):

			if not line:
				continue

			chr, start, gt, svtype, svlen, end, sources, nsvt, genes, ann, svmax, svsum, svtop5, svtop10, svmean, dvg = line

			key = (chr, start, end)
			vdict[key] = (svlen, svtype)

	return vdict

def make_header(samples):
	fields = ["CHR", "START", "END", "N_SAMPLES", "LIST", "GENES", "LONGEST_SVTYPE", "SVLEN", "SVSCORE_MAX", "SVSCORE_SUM", "SVSCORE_TOP5", "SVSCORE_TOP10", "SVSCORE_MEAN", "DGV", "EXONS_SPANNED"]
	header = ",".join(fields)

	for s in samples:
		header = header + "," + s

	for s in samples:
		header = header + "," + "{}_details".format(s)

	return header + "\n"

def make_sample_list_index(samples, interval):

	'''
		LIST column
		e.g. 1;2;3;4;5;6;7
	'''

	index = []
	isthere_array = []

	for i, val in enumerate(samples):
		if val in interval:
			index.append(str(i+1))
			isthere_array.append(str(1))
		else:
			isthere_array.append(str(0))

	#print(index)
	return (";".join(index), isthere_array)

def get_longest_svtype(samples):

	'''
		Determine the longest structural variant, then return its annotation
	'''

	longest = -1	#will always be overwritten on first loop iteration since we are using abs()
	longest_svtype = ""

	for s in samples:
		for variant in samples[s]:

			chr, start, end = variant[0]
			svlen, svtype = variant[1]

			length = abs(int(svlen))

			if length > longest:
				longest = length
				longest_svtype = svtype

	return longest_svtype

def make_sample_details(samples, interval):

	'''
		SAMPLENAME_details column
		e.g. 1:10334731-10334817:DEL
	'''

	all_samp_details = []

	for s in samples:

		samp_details = []

		if s in interval:
			for variant in interval[s]:
				chr, start, end = variant[0]
				svtype = variant[1][1]

				samp_details.append('{}:{}-{}:{}'.format(chr, start, end, svtype))
			all_samp_details.append(';'.join(samp_details))
		else:
			all_samp_details.append("NA")

	return all_samp_details

def write_results(samples, all_sv, all_sv_info, overlapping_exon_count, outfile_name):

	'''
		A lot of string manipulation to generate the final CSV file line by line
	'''

	with open(outfile_name, "w") as out:

		out.write(make_header(samples))

		for key in sorted(all_sv.iterkeys()):

			chr, start, end = key

			index, isthere = make_sample_list_index(samples, all_sv[key])
			nsamples = str(len(all_sv[key]))
			svtype = get_longest_svtype(all_sv[key])
			samp_details = make_sample_details(samples, all_sv[key])

			svlen, svmax, svsum, svtop5, svtop10, svmean, dvg, gene = all_sv_info[key]
			n_exon_spanned = str(overlapping_exon_count[key])

			out_line = '{},{},{}\n'.format(",".join([str(chr), str(start), str(end), nsamples, index, gene, svtype, svlen, svmax, svsum, svtop5, svtop10, svmean, dvg, n_exon_spanned]), ",".join(isthere), ",".join(samp_details))
			out.write(out_line)

def csv2bed(input_files):

	# taken from sergey's previous script - crg.sv.merge_family.sh
	# converts each sample csv to a bed file
	# returns list of sample file names with extension removed

	sample_list = []

	for s in input_files:

		if s.endswith('.sv.csv'):
			s=s[:-7]

		sample_list.append(s)
		cmd = "cat {}.sv.csv | sed 1d | awk -F \'\",\"\' -v smpl={} \'{{print $1\"\\t\"$2\"\\t\"$6\"\\t\"smpl}}\' | sed s/\"\\\"\"// | sort -k1,1 -k2,2n > {}.bed".format(s, s, s)
		os.popen(cmd)

	return sample_list

def count_overlapping_exons(sbed, exon_bed):

	'''
		By referencing a list of fixed exon locations in a bed file, determine the number of exons a SV spans
		Returns a dict where each interval maps to the exon count
	'''

	tmp_bed_name = "tmp_interval.bed"

	exon_ref = BedTool(exon_bed)
	sample = BedTool(sbed)

	overlapping_exon_count = {}

	for interval in sample:

		chr, start, end, gene = interval

		# create a temp bed file with 1 line - the interval of interest
		with open(tmp_bed_name, "w") as f:	# overwrites file with each loop
			f.write('{}\t{}\t{}\n'.format(interval.chrom, interval.start, interval.end))

		tmp_bed = BedTool(tmp_bed_name)

		# count overlapping regions
		n_exons = tmp_bed.intersect(exon_ref).count()

		overlapping_exon_count[(chr, start, end)] = n_exons

	os.remove(tmp_bed_name)	#delete temp file

	return overlapping_exon_count

def make_bed_file(interval, bed_dir):
	with open(bed_dir, "w") as f:
		for i in interval:
			f.write(str(i))

def main(exon_bed, outfile_name, input_files):

	all_sv = {}
	all_sv_info = {}
	overlapping_exon_count = {}
	sample_list = csv2bed(input_files)

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
				pass

		a_bed = current_dir + sample_list[i] + ".bed"
		a_csv = sample_list[i] + ".sv.csv"

		save_sv_column_data(a_csv, all_sv_info)	# store SV info from sample being compared

		for j in range(i, len(sample_list)):

			b_bed = current_dir + sample_list[j] + ".bed"
			b_csv = sample_list[j] + ".sv.csv"

			# print("pass {}: {} {}".format(str(i), a_bed, b_bed) )

			b_svtype_svlen_dict = make_svtype_svlen_dict(b_csv)
			leftover_sv = determine_intersect(a_bed, b_bed, b_svtype_svlen_dict, all_sv)

			# store all leftover_sv in tmp dir for processing in next pass
			if i != j: #if A and B are the same bed file, no need to create tmp file
				new_bed = next_dir + sample_list[j] + ".bed"
				make_bed_file(leftover_sv, new_bed)

	# generate temporary bed file containing all intervals
	with open("all_sv.bed", "w") as f:
		for sv in all_sv:
			f.write('{}\t{}\t{}\t.\n'.format(sv[0], sv[1], sv[2]))

	exons_spanned = count_overlapping_exons("all_sv.bed", exon_bed)
	write_results(sample_list, all_sv, all_sv_info, exons_spanned, outfile_name)

if __name__ == '__main__':

	'''
		How to use: Navigate to a folder containing multiple sv.csv files and launch the script.
		Program generates a CSV file named after outfile variable.

		python crg.intersect_sv_reports.py -exon_bed="protein_coding_genes.exons.bed" -o="180.sv.family.csv" -i 180_123.sv.csv 180_444.sv.csv
	'''

	parser = argparse.ArgumentParser(description='Generates a clinical report in CSV format for structural variants')
	parser.add_argument('-exon_bed', default="/home/naumenko/Desktop/reference_tables/protein_coding_genes.exons.fixed.sorted.bed", help='BED file containing fixed exon positions', required=True)
	parser.add_argument('-o', help='Output file name e.g. 200.sv.family.csv', required=True, type=str)
	parser.add_argument('-i', nargs='+', help='Input file names including .sv.csv extension', required=True)
	args = parser.parse_args()

	print('crg.sv.merge_family.py started processing on ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))
	main(args.exon_bed, args.o, args.i)
	print('crg.sv.merge_family.py finished on ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))
