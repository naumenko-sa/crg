# crg
clinical reseach genome scripts

# Steps:
1. Create project(=case=family) dir:\
`mkdir -p project/input`
2. Copy/symlink input files to project/input: 
- project_sample.bam
- project_sample_1.fq.gz
- project_sample_2.fq.gz
3. `crg.prepare_bcbio_run.sh project`
4. `qsub ~/cre/bcbio.pbs -v project=project`\
or for multiple projects create list of projects in projects.txt and run\
`qsub -t 1-N ~/cre/bcbio.array.pbs`\
where N = number of projects in the current dir.
5. `qsub ~/cre/cre.sh -v family=<project>,cleanup=1,make_report=0,type=wgs`

## AnnotSV
[AnnotSV](http://lbgi.fr/AnnotSV/) must be set up as apart of the local environment to generate family level reports. Users should set FeaturesOverlap and SVtoAnnOverlap to 50 in the configFile. Because these scripts group SV's which have a 50% recipricol overlap, annotation should follow a similar rule.

## Individual sample report columns:
=======
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
crg.intersect_sv_reports.py is a script which groups then annotates structural variants across multiple samples. 

Grouping is useful when analyzing families as most structural variants should be similar and conserved across samples. The criteria for grouping is defined as a minimum of a 50% recipricol overlap with an arbitrary "reference" interval. This criteria gaurentees that grouped structural variants are of similar size and position.

Annotation is largely accomplished by using a tool called [AnnotSV](http://lbgi.fr/AnnotSV/) made by [Véronique Geoffroy](https://www.researchgate.net/profile/Veronique_Geoffroy2). DGV, DDD columns come directly from AnnotSV.
Information on each of these columns can be found on the [AnnotSV website](http://lbgi.fr/AnnotSV/annotations).

The script produces a CSV file which can be analyzed using spreadsheet software.

## Family report columns:
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
- SAMPLENAME: does this sample have an overlapping SV in it? (0,1)
- SAMPLENAME_details: what are the SV's in this sample which overlap with the reference?

## To generate a family level report:
```python crg.intersect_sv_reports.py -exon_bed=/path/to/protein_coding_genes.exons.fixed.sorted.bed -o=output_family_report_name.csv -i sample1.sv.csv sample2.sv.csv sample3.sv.csv```
