import os
from pybedtools import BedTool

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