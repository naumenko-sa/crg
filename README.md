# crg
clinical reseach genome scripts

# Individual sample report columns:
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
- SVscores: [SVscore-github](https://github.com/lganel/SVScore),[SVscore-article](https://academic.oup.com/bioinformatics/article/33/7/1083/2748212)
  - SVSCOREMAX
  - SVSCORESUM
  - SVSCORETOP5
  - SVSCORETOP10
  - SVSCOREMEAN
- DGV: frequency in DGV, 8000 WGS

# Family report

crg.intersect_sv_reports.py is a script which groups structural variants across multiple samples. This is especially useful when analyzing members
of the same family as many mutations should be similar and conserved. The criteria for grouping is a minimum of a 50% recipricol overlap with an
arbitrary "reference" structural variant's interval. The script produces a CSV file which can be analyzed using spreadsheet software to prioritize specific structural variants.

# Family report columns:
Includes all of the columns above, except SOURCES, NUM_SVTOOLS, SVTYPE and ANN, in addition to:
- N_SAMPLES: number of samples a reference interval overlaps with
- LONGEST_SVTYPE: SVTYPE taken from the longest overlapping SV
- EXONS_SPANNED: number of exons a SV affects
- DECIPHER_LINK
- SAMPLENAME: does this sample have an overlapping SV in it? (0,1)
- SAMPLENAME_details: what are the SV's in this sample which overlap with the reference?

## To generate a family level report:
```python crg.intersect_sv_reports.py -exon_bed=/path/to/protein_coding_genes.exons.fixed.sorted.bed -o=output_family_report_name.csv -i sample1.sv.csv sample2.sv.csv sample3.sv.csv```
