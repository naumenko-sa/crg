# 1. Installation
1. Clone [crt](https://github.com/naumenko-sa/crt), [cre](https://github.com/naumenko-sa/cre), and [bioscripts](https://github.com/naumenko-sa/bioscripts) to
`~/crt`, `~/cre`, and `~/bioscripts`.
2. Set PATH: `export PATH=~/crt/scripts:~/cre/scripts:~/bioscripts/scripts` in ~/.bash_profile.
3. Install [bcbio_nextgen](https://github.com/bcbio/bcbio-nextgen), set PATH and PYTHONPATH:
```
export PATH=/path/bcbio/anaconda/bin:$PATH
export PYTHONPATH=/path/bcbio/anaconda/lib/python2.7
```

# 2. Create a bed file for prioritization of small and structural variants
1. Good sources of gene panels are: [Panel App from Genomics England](https://panelapp.genomicsengland.co.uk/) and [Gene panel app from iobio.io](https://genepanel.iobio.io/).
2. If gene list comes from Phenotips:\
```Rscript ~/bioscripts/genes.R phenotips_hpo2gene_coordinates phenotips_hpo.tsv```. Stringr should be >=1.4.
3. Or use [genes.R](https://github.com/naumenko-sa/bioscripts/blob/master/genes.R) in a custom case
4. Some genes might be missing (they don't have ENS IDs in phenotips tsv file, they are reported by script, you can try ~/cre/data/missing_genes_grch37.bed or GeneCards/Ensembl resources to find them).
5. Sort and merge with bedtools\
```
	bedtools sort -i unsorted.bed > sorted.bed
	bedtools merge -i sorted.bed > project.bed
```
6. result is project.bed

# 3. Align reads vs GRCh37 reference with decoy
1. Create a project(=case=family) dir:\
`mkdir -p project/input`
2. Copy/symlink input file(s) to project/input: project_sample.bam, or project_sample_1.fq.gz and project_sample_2.fq.gz
3. Create bcbio project:\
`crg.prepare_bcbio_run.sh project align_decoy`
4. Run bcbio project:\
`qsub ~/cre/scripts/bcbio.pbs -v project=project`\
5. For multiple projects create a list of projects in projects.txt and run:\
`qsub -t 1-N ~/cre/scripts/bcbio.array.pbs`\
where N = number of projects.
6. To speed up the process, run one project per sample.

# 4. Remove decoy reads
`qsub ~/cre/cre.bam.remove_decoy_reads.sh -v bam=$bam`.\
Keep original bam with decoy reads to store all data.\
Some SV callers (manta) are sensitive to reads mapped to decoy even with one mate.

# 5. Call small variants
1. Create a project dir:\
`mkdir -p project/input`
2. Symlink bam file(s) from step1 to project/input: project_sample.bam Small variant calling is not sensitive to the presense of decoy reads.
3. Create bcbio project:\
`crg.prepare_bcbio_run.sh project small_variants`
4. Run bcbio:\
`qsub ~/cre/bcbio.pbs -v project=project`
5. Clean up bcbio project:\
`qsub ~/cre/cre.sh -v family=<project>,cleanup=1,make_report=0,type=wgs`

# 6. Create excel reports for small variants.
1. coding report:\
`qsub ~/cre/cre.sh -v family=project`
2. noncoding variants for gene panels:\
	2.1 subset variants:\
	`bedtools intersect --header -a project-ensemble.vcf.gz -b panel.bed > project.panel.vcf.gz`\
	2.2 reannotate variants in panels and create gemini.db:\
	`qsub ~/crg/crg.vcf2cre.sh -v original_vcf=project.panel.vcf.gz,project=project`\
	**note crg.vcf2cre.sh**  not cre.vcf2cre.sh: annotations for WGS are different  
	2.3 build report:\
	`qsub ~/cre/cre.sh -v family=project,type=wgs`
3. noncoding variants for gene panels with flank\
	3.1 modify bed file, add 100k bp to each gene start and end:\
	`cat panel.bed | awk -F "\t" '{print $1"\t"$2-100000"\t"$3+100000}'`\
	3.2 proceed as for noncoding small variant report
4. de-novo variants for trios

# 7. Call structural variants (in parallel with step 3)
1. MetaSV uses spades (a genome assembler) for every SV, making bcbio run computationally intensive. To speed up call samples individually.\
sv_regions parameter in bcbio works only for wham.
2. Create project dir:\
`mkdir -p project/input`
3. Symlink a bam file. from step 2 to project/input: project_sample.bam.
4. Create bcbio project:\
`crg.prepare_bcbio_run.sh project sv`\
5. Run bcbio:\
`qsub ~/cre/bcbio.pbs -v project=project`

# 8. Excel report for structural variants 
1. [Report columns](https://docs.google.com/document/d/1o870tr0rcshoae_VkG1ZOoWNSAmorCZlhHDpZuZogYE/edit?usp=sharing)
2. `cd project/sv`
3. Generate per sample reports: \
3.1 `cd <sample_dir>`\
3.2 `qsub ~/crg/crg.sv.prioritize.sh -v case=<project>,panel=<panel.bed>`\
Also outputs sample.tsv file for annotation in TCAG with DGV frequencies.\
3.3 (Optional) `crg.sv.prioritize.sh <project> panel.bed tcag_annotated_file.tsv` to incorporate DGV frequencies.
4. Combine reports from multiple samples:\
4.1 Gather each report from the previous step in to a single directory.\
4.3 `crg.intersect_sv_reports.sh project`.

## AnnotSV
[AnnotSV](http://lbgi.fr/AnnotSV/) must be set up as apart of the local environment to generate family level reports. Users should set FeaturesOverlap and SVtoAnnOverlap to 50 in the configFile. Because these scripts group SV's which have a 50% recipricol overlap, annotation should follow a similar rule.

# Report columns:
- CHR
- POS
- GT
- SVTYPE
- SVLEN
- END
- SOURCES: which programs called the event
- NUM_SVTOOLS: how many programs supported the event
- GENES: genes overlapping the event (mostly one gene)
- ANN: raw annotation from VEP
- SVscores: [SVscore-github](https://github.com/lganel/SVScore), [SVscore-article](https://academic.oup.com/bioinformatics/article/33/7/1083/2748212)
  - SVSCOREMAX
  - SVSCORESUM
  - SVSCORETOP5
  - SVSCORETOP10
  - SVSCOREMEAN
- DGV: frequency in DGV, 8000 WGS

# Family report
crg.intersect_sv_reports.py generates a report summarizing structural variants across several samples. It groups structural variants of similar size and position in to a single "reference" structural variant. Grouping is useful when analyzing families as most structural variants should be similar and conserved across samples.

The script produces a CSV file which can be analyzed using spreadsheet software.

## Family report requirements:

### ~/gene_data directory containing the following files:
	HGMD=${HOME}/gene_data/HGMD_2018/hgmd_pro.db
	EXON_BED=${HOME}/gene_data/protein_coding_genes.exons.fixed.sorted.bed
	HPO=${HOME}/gene_data/HPO_2018/${FAMILY_ID}_HPO.txt
	EXAC=${HOME}/gene_data/ExAC/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt
	OMIM=${HOME}/gene_data/OMIM_2018-11-01/genemap2.txt

### AnnotSV
[AnnotSV](http://lbgi.fr/AnnotSV/) must be set up as apart of the local environment to generate a family report. Users should set FeaturesOverlap and SVtoAnnOverlap to 50 in the configFile. Because these scripts group SV's which have a 50% recipricol overlap, annotation should follow a similar rule.

DGV and DDD columns are annotated by [AnnotSV](http://lbgi.fr/AnnotSV/).

## Family report columns:
[In-depth column descriptions](https://docs.google.com/document/d/1o870tr0rcshoae_VkG1ZOoWNSAmorCZlhHDpZuZogYE/edit#)

Includes all of the columns above, except SOURCES, NUM_SVTOOLS, SVTYPE and ANN, in addition to:
- N_SAMPLES: number of samples a reference interval overlaps with
- LONGEST_SVTYPE: SVTYPE taken from the longest overlapping SV
- EXONS_SPANNED: number of exons a SV affects
- DECIPHER_LINK: hyperlink to DECIPHER website for the reference interval
- DGV_GAIN_IDs
- DGV_GAIN_n_samples_with_SV
- DGV_GAIN_n_samples_tested
- DGV_GAIN_Frequency
- DGV_LOSS_IDs
- DGV_LOSS_n_samples_with_SV
- DGV_LOSS_n_samples_tested
- DGV_LOSS_Frequency
- DDD_SV
- DDD_DUP_n_samples_with_SV
- DDD_DUP_Frequency
- DDD_DEL_n_samples_with_SV
- DDD_DEL_Frequency
- OMIM MIM# 
- OMIM INHERITANCE
- OMIM DESCRIPTION
- SYNZ
- MISZ
- PLI
- GENES_IN_HGMD
- HGMD_DISEASE
- HGMD_TAG
- HGMD_DESCRIPTION
- HGMD_COMMENT
- HGMD_JOURNAL_INFO
- HGMD_GROSS_INSERTION: gross (>20bp) insertion events in this gene that have been observed in HGMD
- HGMD_GROSS_DUPLICATION: gross duplication events in this gene that have been observed in HGMD
- HGMD_GROSS_DELETION: gross deletion events in this gene that have been observed in HGMD
- HGMD_COMPLEX_VARIATION: complex variations (combination of indels, translocations, SNP, fusions, inversions) in this gene that have been observed in HGMD
- SAMPLE: does this sample have an overlapping SV in it? (0,1)
- SAMPLE_details: what are the SV's in this sample which overlap with the reference?
- SAMPLE_GENOTYPE

## Result dir structure:
project(family)_ID:
- bcbio-align: config and final dirs from bcbio align-decoy run
- bcbio-small-variants: bcbio configs and vcfs from bcbio small variant, output of cre for coding report
- bcbio-sv: SV output from bcbio
- genes: HPO, gene list, bed file
- panel: non-coding report for gene panel cre dir
- panel-flank100k: non-coding report for gene panel +100k flank cre dir
- reports: csv report we send
- tcag: tcag analysis
- bam and bai files (without ready) - in the top directory to for easy access, bams from align-decoy step!

## Use case: compared SV calls from TCAG (ERDS) to MetaSV
```
bcftools view -i 'ALT="<DEL>"' 159_CH0315.pass.vcf.gz | bcftools query -f '%CHROM\t%POS\t%INFO/END\n' -o 159.metasv.del.bed
cat 159_CH0315.erds+_db20171204_20180815_3213_annotated.tsv |  awk '$5 ~/DEL/{print $2"\t"$3"\t"$4}'  > 159.tcag.del.bed
bedtools intersect -a 159.tcag.dup.bed -b 159.metasv.dup.bed -f 0.5 -wo -r | wc -l
```
