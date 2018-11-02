import os
import subprocess
import sqlite3
from pybedtools import BedTool

class GeneAnnotations:
	def __init__(self):
		self.gene = ""

		#DDD fields
		self.DDD_status = ""
		self.DDD_mode = ""
		self.DDD_consequence = ""
		self.DDD_disease = ""
		self.DDD_pmids = ""

		#HGMD GROSS mutation fields
		self.is_in_hgmd = ""
		self.hgmd_disease = ""
		self.hgmd_tag = ""
		self.hgmd_description = ""
		self.hgmd_comments = ""
		self.hgmd_journal = ""
		self.hgmd_author = ""
		self.hgmd_year = ""
		self.hgmd_pmid = ""

		#OMIM fields
		self.mim_num = ""
		self.mim_inheritance = ""
		self.mim_description = ""

		#pLI scores
		self.synz = ""
		self.misz = ""
		self.pli = ""

class StructuralVariant:
	'''Simple class to hold column values'''

	def __init__(self):
		#Standard BED fields
		self.chr = ""
		self.start = ""
		self.end = ""
		self.gene = ""
		self.genotype = ""

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

		#AnnotSV - DGV Fields
		self.dgv_gain_id = ""
		self.dgv_gain_n_samples_with_sv = ""
		self.dgv_gain_n_samples_tested = ""
		self.dgv_gain_freq = ""

		self.dgv_loss_id = ""
		self.dgv_loss_n_samples_with_sv = ""
		self.dgv_loss_n_samples_tested = ""
		self.dgv_loss_freq = ""

		#AnnotSV - DDD fields
		self.ddd_sv = ""
		self.ddd_dup_n_samples_with_sv = ""
		self.ddd_dup_freq = ""
		self.ddd_del_n_samples_with_sv = ""
		self.ddd_del_freq = ""

		#AnnotSV - OMIM fields
		self.omim = []

		#AnnotSV - pLi scores
		self.synz = ""
		self.misz = ""
		self.pli = ""

		#HGMD SVs observed in the same gene
		self.hgmd_gross_deletion = []
		self.hgmd_gross_duplication = []
		self.hgmd_gross_insertion = []
		self.is_gene_in_hgmd = {}

	def make_decipher_link(self):
		return '=hyperlink("https://decipher.sanger.ac.uk/browser#q/{}:{}-{}")'.format(self.chr, self.start, self.end)

	def key(self):
		return (self.chr, self.start, self.end)

	def interval(self):
		return '{}:{}-{}:{}'.format(self.chr, self.start, self.end, self.svtype)

	def make_omim_column(self):
		return '; '.join(self.omim)

	def add_omim(self, gene, mim_number, inheritance, phenotype):
		self.omim.append('{GENE: %s MIM_NUMBER: %s INHERITANCE: %s PHENOTYPES: %s}' % (gene, mim_number, inheritance, phenotype))

	def make_hgmd_column(self, hgmd_sv):
		return '; '.join('{%s}' % study for study in hgmd_sv)

class StructuralVariantRecords:
	'''
		Holds groupings of StructuralVariant.

		Groupings in the dict grouped_sv follows this structure:
			grouped_sv[ref_interval][samp_name] = [StructuralVariant]

		Annotated information on the ref_interval used to group sample intervals is stored in:
			all_ref_interval_data[(chr, start, end)] = StructuralVariant
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
		if ref_interval not in self.grouped_sv.keys():
			self.grouped_sv[ref_interval] = {}
			self.all_ref_interval_data[ref_interval] = column_data[ref_interval]

		if samp_name not in self.grouped_sv[ref_interval]:
			self.grouped_sv[ref_interval][samp_name] = []

		self.grouped_sv[ref_interval][samp_name].append(column_data[new_interval])

	def make_header(self):
		fields = ["CHR", "START", "END", "N_SAMPLES", "LIST", "GENES", "LONGEST_SVTYPE", "SVLEN", "SVSCORE_MAX", \
		"SVSCORE_SUM", "SVSCORE_TOP5", "SVSCORE_TOP10", "SVSCORE_MEAN", "DGV", "EXONS_SPANNED", "DECIPHER_LINK", \
		"DGV_GAIN_IDs", "DGV_GAIN_n_samples_with_SV", "DGV_GAIN_n_samples_tested", "DGV_GAIN_Frequency", "DGV_LOSS_IDs", \
		"DGV_LOSS_n_samples_with_SV", "DGV_LOSS_n_samples_tested", "DGV_LOSS_Frequency", "DDD_SV", "DDD_DUP_n_samples_with_SV", \
		"DDD_DUP_Frequency", "DDD_DEL_n_samples_with_SV", "DDD_DEL_Frequency", "OMIM {GENE MIM# INHERITANCE DESCRIPTION};", \
		"synZ", "misZ", "pLI", "GENE_IN_HGMD", "HGMD_GROSS_INSERTION", "HGMD_GROSS_DUPLICATION", "HGMD_GROSS_DELETION"]
		fields.extend(self.sample_list)
		fields.extend(["%s_details" % s for s in self.sample_list])
		return "%s\n" % ",".join(fields)

	def make_bed(self, bed_name):
		'''
			Creates a bed file containing all reference intervals
		'''
		with open(bed_name, "w") as f:
			for sv in self.grouped_sv:
				f.write('{}\t{}\t{}\t{}\n'.format(sv[0], sv[1], sv[2], bed_name))

	def write_results(self, outfile_name):
		'''
			A lot of string manipulation to generate the final CSV file line by line
		'''

		def make_sample_list_index(sample_list, interval):
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

		def make_sample_details(sample_list, interval):
			'''
				Makes list for SAMPLENAME_details column
				e.g. 1:10334731-10334817:DEL;1:10334769-10334833:DUP
			'''
			all_samp_details = []
			for s in sample_list:
				samp_details = []

				if s in interval:
					for variant in interval[s]:
						samp_details.append(variant.interval())
					all_samp_details.append(';'.join(samp_details))
				else:
					all_samp_details.append("NA")

			return all_samp_details

		def make_sample_genotype_details(sample_list, interval):
			return ["NA" if sample not in interval.keys() \
			else "HOM" if "HOM" in [variant.genotype for variant in interval.values()] \
			else "HET" \
			for sample in sample_list]

		def get_longest_svtype(interval):
			'''
				Populates LONGEST_SVTYPE column
			'''
			def svlen(type_len):
				return abs(int(type_len[1]))

			svtype_and_svlen = [(variant.svtype, variant.svlen) for sample in interval.values() for variant in sample]
			svtype_and_svlen.sort(key=svlen)
			return svtype_and_svlen[-1][0] #longest SV is last value in sorted list, hence [-1]

		def hgmd_gene_index(gene_dict):
			return u"; ".join(["{}: {}".format(gene, present) for gene, present in gene_dict.items()])

		with open(outfile_name, "w") as out:
			out.write(self.make_header())

			for key in sorted(self.grouped_sv.keys()):

				chr, start, end = key
				ref = self.all_ref_interval_data[key]

				index, isthere = make_sample_list_index(self.sample_list, self.grouped_sv[key])
				nsamples = str(len(self.grouped_sv[key]))
				svtype = get_longest_svtype(self.grouped_sv[key])
				samp_details = make_sample_details(self.sample_list, self.grouped_sv[key])

				out_line = u'{},{},{}\n'.format(",".join([str(chr), str(start), str(end), nsamples, index, ref.gene, svtype, \
				ref.svlen, ref.svmax, ref.svsum, ref.svtop5, ref.svtop10, ref.svmean, ref.dgv, ref.exons_spanned, ref.make_decipher_link(), \
				ref.dgv_gain_id, ref.dgv_gain_n_samples_with_sv, ref.dgv_gain_n_samples_tested, ref.dgv_gain_freq, ref.dgv_loss_id, \
				ref.dgv_loss_n_samples_with_sv, ref.dgv_loss_n_samples_tested, ref.dgv_loss_freq, ref.ddd_sv, ref.ddd_dup_n_samples_with_sv, \
				ref.ddd_dup_freq, ref.ddd_del_n_samples_with_sv, ref.ddd_del_freq, ref.make_omim_column(), ref.synz, ref.misz, ref.pli,\
				hgmd_gene_index(ref.is_gene_in_hgmd), ref.make_hgmd_column(ref.hgmd_gross_insertion), ref.make_hgmd_column(ref.hgmd_gross_duplication), \
				ref.make_hgmd_column(ref.hgmd_gross_deletion)]), ",".join(isthere), ",".join(samp_details))

				out.write(out_line)

	def annotate(self, exon_bed, hgmd_db):
		def calc_exons_spanned(exon_bed):
			'''
				Populates the field: exons_spanned for all reference intervals
			'''
			tmp_bed_name = "tmp_interval.bed"
			tmp_all_sv_bed_name = "tmp_all_sv.bed"

			exon_ref = BedTool(exon_bed)
			self.make_bed(tmp_all_sv_bed_name)
			sample = BedTool(tmp_all_sv_bed_name)

			for interval in sample:

				chr, start, end, gene = interval

				# create a temp bed file with 1 line - the interval of interest
				with open(tmp_bed_name, "w") as f:
					f.write('{}\t{}\t{}\n'.format(chr, start, end))

				tmp_bed = BedTool(tmp_bed_name)
				self.all_ref_interval_data[(chr, start, end)].exons_spanned = str(tmp_bed.intersect(exon_ref).count())

			os.remove(tmp_bed_name)
			os.remove(tmp_all_sv_bed_name)

		def annotsv():
			'''
			Handles DGV, DDD and OMIM annotations
			'''
			all_sv_bed_name = "all_sv.bed"
			annotated = "./{}.annotated.tsv".format(all_sv_bed_name)
			self.make_bed(all_sv_bed_name)

			subprocess.call("$ANNOTSV/bin/AnnotSV -SVinputFile {} -SVinputInfo 1 -outputFile {}".format(all_sv_bed_name, annotated), shell=True)

			with open(annotated, "r") as f:

				next(f) #skip header

				for line in f:
					field = line.rstrip('\n').replace(',', ';').split('\t')
					sv = self.all_ref_interval_data[(field[0], field[1], field[2])]
					if field[4] == "full":
						sv.dgv_gain_id = field[13]
						sv.dgv_gain_n_samples_with_sv = field[14]
						sv.dgv_gain_n_samples_tested = field[15]
						sv.dgv_gain_freq = field[16]

						sv.dgv_loss_id = field[17]
						sv.dgv_loss_n_samples_with_sv = field[18]
						sv.dgv_loss_n_samples_tested = field[19]
						sv.dgv_loss_freq = field[20]

						sv.ddd_sv = field[21]
						sv.ddd_dup_n_samples_with_sv = field[22]
						sv.ddd_dup_freq = field[23]
						sv.ddd_del_n_samples_with_sv =field[24]
						sv.ddd_del_freq = field[25]
					elif field[4] == "split":
						# gene info for omim annotations
						sv.add_omim(field[5], field[38], field[40], field[39]) # gene, mim_number, phenotypes, inheritance

						# pli scores, identical for each sv interval
						sv.synz = field[34]
						sv.misz = field[35]
						sv.pli = field[36]

			os.remove(all_sv_bed_name)
			os.remove(annotated)

		def hgmd(db_path):
			def decode(rows):
				return [ u'|'.join(["NA" if field is None else str(field) if isinstance(field,int) else field.replace(',', '') for field in entry ]) for entry in rows ]

			conn = sqlite3.connect(db_path)
			cur = conn.cursor()

			for sv in self.all_ref_interval_data.values():
				svtype = sv.svtype

				for gene in set(sv.gene.replace('&', ';').split(';')): #set ensures uniqueness - the same gene don't get querried twice
					if gene.isspace() or len(gene) == 0:
						continue

					# Is this gene in the HGMD database?
					cur.execute('SELECT GENE FROM ALLGENES WHERE GENE=?', (gene, ))
					sv.is_gene_in_hgmd[gene] = "0" if cur.fetchone() is None else "1"

					# Look for solved cases in this gene involving structural variants 
					if svtype == "DEL":
						cur.execute('SELECT GENE, DISEASE, TAG, DESCR, COMMENTS, JOURNAL, AUTHOR, YEAR, PMID FROM GROSDEL WHERE GENE=?', (gene, ) )
						sv.hgmd_gross_deletion.extend(decode(cur.fetchall()))
					elif svtype == "INS":
						cur.execute('SELECT GENE, DISEASE, TAG, DESCR, COMMENTS, JOURNAL, AUTHOR, YEAR, PMID FROM GROSINS WHERE GENE=? AND TYPE=?', (gene, 'I'))
						sv.hgmd_gross_insertion.extend(decode(cur.fetchall()))
					elif svtype == "DUP":
						cur.execute('SELECT GENE, DISEASE, TAG, DESCR, COMMENTS, JOURNAL, AUTHOR, YEAR, PMID FROM GROSINS WHERE GENE=? AND TYPE=?', (gene, 'D'))
						sv.hgmd_gross_duplication.extend(decode(cur.fetchall()))

			conn.close()

		calc_exons_spanned(exon_bed)
		annotsv()
		hgmd(hgmd_db)