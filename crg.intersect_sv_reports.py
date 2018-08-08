import os
import csv
import re
import sys

outfile = "overlapping_family_sv.csv"
csv.field_size_limit(sys.maxsize)

def determine_intersect(first_sample_bed, other_sample_bed, other_sample_sv, first_sample_sv):

	# bedtools doc: https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html?highlight=intersect
	# find all SVs in other_sample_bed which overlap a SV in first_sample_bed by >=50% and are overlapped by >=50% by the same SV in first_sample_bed
	out = os.popen("bedtools intersect -a {} -b {} -wa -wb -F 0.5 -f 0.5".format(first_sample_bed, other_sample_bed)).read()

	for l in out.split("\n"):
		if not l:
			continue

		chr, start, end, ref_name, samp_chr, samp_start, samp_end, samp_name = l.split("\t")

		interval = (chr, start, end)
		other_interval = (samp_chr, samp_start, samp_end)

		if samp_name not in first_sample_sv[interval]:
			first_sample_sv[interval][samp_name] = []

		first_sample_sv[interval][samp_name].append((other_interval, (other_sample_sv[other_interval])))

def parse_first_sample_csv(scsv):
	with open(scsv) as f:

		next(f)	#skip the header
		vdict = {}
		variant_info = {}

		for line in csv.reader(f, delimiter=",", quotechar="\""):

			if not line:
				continue

			chr, start, gt, svtype, svlen, end, sources, nsvt, genes, ann, svmax, svsum, svtop5, svtop10, svmean, dvg = line
			key = (chr, start, end)

			vdict[key] = {}
			variant_info[key] = (svlen, svmax, svsum, svtop5, svtop10, svmean, dvg)

	return vdict, variant_info

def parse_sample_csv(scsv):

	#read and create svtype dict from sample csv
	#vdict format: vdict[(chr, start, end)] = svtype

	with open(scsv) as f:

		next(f)	#skip the header
		vdict = {}

		for line in csv.reader(f, delimiter=",", quotechar="\""):

			if not line:
				continue

			chr, start, gt, svtype, svlen, end, sources, nsvt, genes, ann, svmax, svsum, svtop5, svtop10, svmean, dvg = line

			key = (chr, start, end)
			vdict[key] = svtype

	return vdict

def make_header(samples):
	fields = ["CHR","START", "END", "N_SAMPLES", "LIST", "LONGEST_SVTYPE", "SVLEN", "SVSCORE_MAX", "SVSCORE_SUM", "SVSCORE_TOP5", "SVSCORE_TOP10", "SVSCORE_MEAN", "DVG"]
	header = ",".join(fields)

	for s in samples:
		header = header + "," + s

	for s in samples:
		header = header + "," + "{}_details".format(s)

	return header + "\n"

def make_sample_list_index(samples, interval):

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
	#determine the longest structural variant
	#then return its annotation

	longest = -1	#will always be overwritten on first loop iteration since we are using abs()
	svtype = ""

	for s in samples:
		for variant in samples[s]:
			chr, start, end = variant[0]

			length = abs(int(end) - int(start))

			if length > longest:
				length = longest
				svtype = variant[1]
				#print("variant: " + variant[1])

	return svtype

def make_sample_details(samples, interval):

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

def write_results(samples, first_sample_sv, first_sample_sv_info):

	with open(outfile, "w") as out:

		out.write(make_header(samples))

		for key in sorted(first_sample_sv.iterkeys()):

			chr, start, end = key

			index, isthere = make_sample_list_index(samples, first_sample_sv[key])
			nsamples = str(len(first_sample_sv[key]))
			svtype = get_longest_svtype(first_sample_sv[key])
			samp_details = make_sample_details(samples, first_sample_sv[key])

			svlen, svmax, svsum, svtop5, svtop10, svmean, dvg = first_sample_sv_info[key]

			out_line = '{},{},{}\n'.format(",".join([str(chr), str(start), str(end), nsamples, index, svtype, svlen, svmax, svsum, svtop5, svtop10, svmean, dvg]), ",".join(isthere), ",".join(samp_details))
			out.write(out_line)

def csv2bed():
	#taken from sergey's previous script - crg.sv.merge_family.sh
	'''ls *.sv.csv | sed s/.sv.csv// > samples.txt

		for sample in `cat samples.txt`
		do
		    cat $sample.sv.csv | sed 1d | awk -F '","' -v smpl=$sample '{print $1"\t"$2"\t"$6"\t"smpl}' | sed s/"\""// | sort -k1,1 -k2,2n > $sample.bed
		done'''

	#creates sample.txt
	#converts each sample csv to a bed file
	os.popen("ls *.sv.csv | sed s/.sv.csv// > samples.txt")

	with open("samples.txt") as f:
		for s in f:
			s = s.strip()
			cmd = "cat {}.sv.csv | sed 1d | awk -F \'\",\"\' -v smpl={} \'{{print $1\"\\t\"$2\"\\t\"$6\"\\t\"smpl}}\' | sed s/\"\\\"\"// | sort -k1,1 -k2,2n > {}.bed".format(s, s, s)
			os.popen(cmd)

def main():
	csv2bed()

	sample_list = []
	first_sample_sv = {}
	first_sample_bed = ""
	first_sample_sv_info = {}

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

	write_results(sample_list, first_sample_sv, first_sample_sv_info)

if __name__ == '__main__':
	# Launch script within a folder containing multiple sv.csv files
	print('crg.sv.merge_family.py started processing.')
	main()
	print('crg.sv.merge_family.py is done processing!')
