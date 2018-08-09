import os
import csv
import re
import sys
import argparse

from datetime import datetime
from pybedtools import BedTool

outfile = "overlapping_family_sv.csv"
csv.field_size_limit(sys.maxsize)

def determine_intersect(first_sample_bed, other_sample_bed, other_sample_sv, first_sample_sv):
	
	'''
		Populate first_sample_sv dict with overlapping SVs from another bed file.
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

		if samp_name not in first_sample_sv[interval]:
			first_sample_sv[interval][samp_name] = []

		first_sample_sv[interval][samp_name].append((other_interval, (other_sample_sv[other_interval])))

def parse_first_sample_csv(scsv):

	'''
		Parse the first sample's CSV and set up data structures to store interval data.

		variant_info holds SV_SCOREs and other column info from the CSV file
		vdict is to be later used to hold variant info from other samples which fall within a specific interval
	'''
	with open(scsv) as f:

		vdict = {}
		variant_info = {}

		next(f)	# skip the header

		for line in csv.reader(f, delimiter=",", quotechar="\""):

			if not line:
				continue

			chr, start, gt, svtype, svlen, end, sources, nsvt, genes, ann, svmax, svsum, svtop5, svtop10, svmean, dvg = line
			key = (chr, start, end)

			vdict[key] = {}
			variant_info[key] = (svlen, svmax, svsum, svtop5, svtop10, svmean, dvg)

	return vdict, variant_info

def parse_sample_csv(scsv):

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
			vdict[key] = svtype

	return vdict

def make_header(samples):
	fields = ["CHR", "START", "END", "N_SAMPLES", "LIST", "LONGEST_SVTYPE", "SVLEN", "SVSCORE_MAX", "SVSCORE_SUM", "SVSCORE_TOP5", "SVSCORE_TOP10", "SVSCORE_MEAN", "DVG", "EXONS_SPANNED"]
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
	svtype = ""

	for s in samples:
		for variant in samples[s]:
			chr, start, end = variant[0]

			length = abs(int(end) - int(start))

			if length > longest:
				length = longest
				svtype = variant[1]

	return svtype

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
				svtype = variant[1]

				samp_details.append('{}:{}-{}:{}'.format(chr, start, end, svtype))
			all_samp_details.append(';'.join(samp_details))
		else:
			all_samp_details.append("NA")

	return all_samp_details

def write_results(samples, first_sample_sv, first_sample_sv_info, overlapping_exon_count):

	'''
		A lot of string manipulation to generate the final CSV file line by line
	'''

	with open(outfile, "w") as out:

		out.write(make_header(samples))

		for key in sorted(first_sample_sv.iterkeys()):

			chr, start, end = key

			index, isthere = make_sample_list_index(samples, first_sample_sv[key])
			nsamples = str(len(first_sample_sv[key]))
			svtype = get_longest_svtype(first_sample_sv[key])
			samp_details = make_sample_details(samples, first_sample_sv[key])

			svlen, svmax, svsum, svtop5, svtop10, svmean, dvg = first_sample_sv_info[key]
			n_exon_spanned = str(overlapping_exon_count[key])

			out_line = '{},{},{}\n'.format(",".join([str(chr), str(start), str(end), nsamples, index, svtype, svlen, svmax, svsum, svtop5, svtop10, svmean, dvg, n_exon_spanned]), ",".join(isthere), ",".join(samp_details))
			out.write(out_line)

def csv2bed():
	# taken from sergey's previous script - crg.sv.merge_family.sh

	'''
		ls *.sv.csv | sed s/.sv.csv// > samples.txt

		for sample in `cat samples.txt`
		do
		    cat $sample.sv.csv | sed 1d | awk -F '","' -v smpl=$sample '{print $1"\t"$2"\t"$6"\t"smpl}' | sed s/"\""// | sort -k1,1 -k2,2n > $sample.bed
		done
	'''

	# creates sample.txt
	# converts each sample csv to a bed file

	os.popen("ls *.sv.csv | sed s/.sv.csv// > samples.txt")

	with open("samples.txt") as f:
		for s in f:
			s = s.strip()
			cmd = "cat {}.sv.csv | sed 1d | awk -F \'\",\"\' -v smpl={} \'{{print $1\"\\t\"$2\"\\t\"$6\"\\t\"smpl}}\' | sed s/\"\\\"\"// | sort -k1,1 -k2,2n > {}.bed".format(s, s, s)
			os.popen(cmd)

def find_overlapping_exons(sv_dict, sbed, exon_bed):

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

def main(exon_bed):

	csv2bed()

	sample_list = []
	first_sample_sv = {}
	first_sample_bed = ""
	first_sample_sv_info = {}
	overlapping_exon_count = {}

	with open("samples.txt") as f:
		for i, s in enumerate(f):

			s = s.strip()	# newline

			sbed = s + ".bed"
			scsv = s + ".sv.csv"
			sample_list.append(s)

			if i == 0:
				first_sample_sv, first_sample_sv_info = parse_first_sample_csv(scsv)
				first_sample_bed = sbed

			other_sample_sv = parse_sample_csv(scsv)
			determine_intersect(first_sample_bed, sbed, other_sample_sv, first_sample_sv)

	overlapping_exon_count = find_overlapping_exons(first_sample_sv, first_sample_bed, exon_bed)

	write_results(sample_list, first_sample_sv, first_sample_sv_info, overlapping_exon_count)

if __name__ == '__main__':

	'''
		How to use: Navigate to a folder containing multiple sv.csv files and launch the script.
		Program generates a CSV file named after outfile variable.

		python crg.intersect_sv_reports.py -exon_bed="/absolute/path/to/protein_coding_genes.exons.bed"
	'''

	parser = argparse.ArgumentParser(description='Generate a clinical report in CSV format for structural variants')
	parser.add_argument('-exon_bed', default="/home/naumenko/Desktop/reference_tables/protein_coding_genes.exons.fixed.sorted.bed", help='BED file containing fixed exon positions')

	args = parser.parse_args()

	print('crg.sv.merge_family.py started processing on ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))
	main(args.exon_bed)
	print('crg.sv.merge_family.py finished on ' + datetime.now().strftime("%Y-%m-%d_%H:%M:%S.%f"))
