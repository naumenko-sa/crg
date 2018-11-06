import os
import subprocess
import sqlite3
from pybedtools import BedTool

class GeneAnnotations:
	def __init__(self, gene_name):
		self.gene_name = gene_name

		#DDD fields
		self.DDD_status = ""
		self.DDD_mode = ""
		self.DDD_consequence = ""
		self.DDD_disease = ""
		self.DDD_pmids = ""

		self.is_in_hgmd = False
		#HGMD GROSS mutation fields
		self.hgmd_disease = []
		self.hgmd_tag = []
		self.hgmd_description = []
		self.hgmd_comments = []
		self.hgmd_journal = []
		self.hgmd_author = []
		self.hgmd_year = []
		self.hgmd_pmid = []

		#OMIM fields
		self.mim_num = ""
		self.mim_inheritance = ""
		self.mim_description = ""

		#pLI scores
		self.synz = ""
		self.misz = ""
		self.pli = ""

	def add_hgmd_anno(self, hgmd_annots):
		for row in hgmd_annots:
			disease, tag, description, comments, journal, author, year, pmid = row
			self.hgmd_disease.append(disease)
			self.hgmd_tag.append(tag)
			self.hgmd_description.append(description)
			self.hgmd_comments.append(comments)
			self.hgmd_journal.append(journal)
			self.hgmd_author.append(author)
			self.hgmd_year.append(year)
			self.hgmd_pmid.append(pmid)

class StructuralVariant:
	def __init__(self, chr, start, end, svtype, genotype, svlen, svsum, svmax, svtop5, svtop10, svmean, dgv):
		self.chr = chr
		self.start = start
		self.end = end
		self.genotype = genotype

		# GeneAnnotation values 
		self.genes = {}

		#SVSCORES fields
		self.svlen = svlen
		self.svtype = svtype
		self.svmax = svmax
		self.svsum = svsum
		self.svtop5 = svtop5
		self.svtop10 = svtop10
		self.svmean = svmean

		# Custom implementation
		self.exons_spanned = ""

		#Annotated by TCAG
		self.dgv = dgv

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

	def make_decipher_link(self):
		return '=hyperlink("https://decipher.sanger.ac.uk/browser#q/{}:{}-{}")'.format(self.chr, self.start, self.end)

	def key(self):
		return (self.chr, self.start, self.end)

	def make_interval_string(self):
		return '{}:{}-{}:{}'.format(self.chr, self.start, self.end, self.svtype)
	
	def add_gene(self, gene_name):
		if gene not in self.genes.keys():
			self.genes[gene_name] = GeneAnnotations(gene_name)
		else:
			raise ValueError('Attempted to add %s twice to %s' % (gene_name, self.make_interval_string()))

	def make_hgmd_gene_index(self):
		return ';'.join([gene.gene_name for gene in self.genes if gene.is_in_hgmd])

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
		"synZ", "misZ", "pLI", "GENES_IN_HGMD", "HGMD_GROSS_INSERTION", "HGMD_GROSS_DUPLICATION", "HGMD_GROSS_DELETION"]
		fields.extend(self.sample_list)
		fields.extend(["%s_SV_details" % s for s in self.sample_list])
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
				Makes strings for LIST and SAMPLENAME column
				e.g. 1;2;3;4;5;6;7 and 1,1,1,1,1,1,1
			'''
			return ";".join([str(i+1) if sample in interval for i, sample in enumerate(self.sample_list)]), \
			",".join(["1" if sample in interval else "0" for sample in sample_list])

		def make_sample_sv_details(sample_list, interval):
			'''
				Makes list for SAMPLENAME_details column
				e.g. 1:10334731-10334817:DEL;1:10334769-10334833:DUP
			'''
			return ','.join(["NA" if sample not in interval \
			else ";".join([variant.make_interval_string() for variant in interval[sample]]) \
			for sample in sample_list ])

		def make_sample_genotype_details(sample_list, interval):
			return ','.join(["NA" if sample not in interval \
			else "HOM" if "HOM" in [variant.genotype for variant in interval.values()] \
			else "HET" \
			for sample in sample_list])

		def get_longest_svtype(interval):
			'''
				Returns SVTYPE of largest SV in a grouping
			'''
			def svlen(type_len):
				return abs(int(type_len[1]))

			svtype_and_svlen = [(variant.svtype, variant.svlen) for sample in interval.values() for variant in sample]
			svtype_and_svlen.sort(key=svlen)
			return svtype_and_svlen[-1][0] #longest SV is last value in sorted list, hence [-1]

		with open(outfile_name, "w") as out:
			out.write(self.make_header())

			for key in sorted(self.grouped_sv.keys()):

				#TODO: Rewrite output file formatting for GeneAnnotations

				chr, start, end = key
				ref = self.all_ref_interval_data[key]

				sample_list_index, sample_list_isthere_index = make_sample_list_index(self.sample_list, self.grouped_sv[key])
				n_samples = str(len(self.grouped_sv[key]))
				svtype = get_longest_svtype(self.grouped_sv[key])
				samp_sv_details = make_sample_sv_details(self.sample_list, self.grouped_sv[key])

				out_line = u'{},{},{}\n'.format(",".join([str(chr), str(start), str(end), n_samples, sample_list_index, ref.gene, svtype, \
				ref.svlen, ref.svmax, ref.svsum, ref.svtop5, ref.svtop10, ref.svmean, ref.dgv, ref.exons_spanned, ref.make_decipher_link(), \
				ref.dgv_gain_id, ref.dgv_gain_n_samples_with_sv, ref.dgv_gain_n_samples_tested, ref.dgv_gain_freq, ref.dgv_loss_id, \
				ref.dgv_loss_n_samples_with_sv, ref.dgv_loss_n_samples_tested, ref.dgv_loss_freq, ref.ddd_sv, ref.ddd_dup_n_samples_with_sv, \
				ref.ddd_dup_freq, ref.ddd_del_n_samples_with_sv, ref.ddd_del_freq, ref.make_omim_column(), ref.synz, ref.misz, ref.pli, \
				ref.make_hgmd_gene_index(), ref.make_hgmd_column(ref.hgmd_gross_insertion), ref.make_hgmd_column(ref.hgmd_gross_duplication), \
				ref.make_hgmd_column(ref.hgmd_gross_deletion)]), \
				sample_list_isthere_index, samp_sv_details)

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
						#TODO: Rewrite for GeneAnnotations

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
				return [ ["NA" if field is None else \
				str(field) if isinstance(field,int) else \
				field.replace(',', ';') \
				for field in row ] \
				for row in rows ]

			conn = sqlite3.connect(db_path)
			cur = conn.cursor()

			for sv in self.all_ref_interval_data.values():
				svtype = sv.svtype

				for gene_anno in sv.genes:
					gene_name = gene_anno.gene_name
					cur.execute('SELECT GENE FROM ALLGENES WHERE GENE=?', (gene_name, ))

					if cur.fetchone() is not None:
						gene_anno.is_in_hgmd = True
						# Look for solved cases in this gene involving structural variants 
						if svtype == "DEL":
							cur.execute('SELECT DISEASE, TAG, DESCR, COMMENTS, JOURNAL, AUTHOR, YEAR, PMID FROM GROSDEL WHERE GENE=?', (gene_name, ) )
							gene_anno.add_hgmd_anno(decode(cur.fetchall()))
						elif svtype == "INS":
							cur.execute('SELECT DISEASE, TAG, DESCR, COMMENTS, JOURNAL, AUTHOR, YEAR, PMID FROM GROSINS WHERE GENE=? AND TYPE=?', (gene_name, 'I'))
							gene_anno.add_hgmd_anno(decode(cur.fetchall()))
						elif svtype == "DUP":
							cur.execute('SELECT DISEASE, TAG, DESCR, COMMENTS, JOURNAL, AUTHOR, YEAR, PMID FROM GROSINS WHERE GENE=? AND TYPE=?', (gene_name, 'D'))
							gene_anno.add_hgmd_anno(decode(cur.fetchall()))
					else:
						gene_anno.is_in_hgmd = False

			conn.close()

		calc_exons_spanned(exon_bed)
		annotsv()
		hgmd(hgmd_db)