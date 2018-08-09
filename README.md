# crg
clinical reseach genome scripts

# Steps:
1. Create project(=case=family) dir: mkdir -p project/input.
2. Copy/symlink input files in project/input: project_sample.bam, project_sample_1.fq.gz,project_sample_2.fq.gz.
3. crg.prepare_bcbio_run.sh project
4. qsub ~/cre/bcbio.pbs -v project=project
or for multiple projects create projects.txt and run

qsub -t 1-N ~/cre/bcbio.array.pbs

where N = number of projects in the current dir.

5. qsub ~/cre/cre.sh -v family=<project>,cleanup=1,make_report=0,type=wgs - cleaning up temporary files

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
- SVscores: [SVscore-github](https://github.com/lganel/SVScore),[SVscore-article](https://academic.oup.com/bioinformatics/article/33/7/1083/2748212)
  - SVSCOREMAX
  - SVSCORESUM
  - SVSCORETOP5
  - SVSCORETOP10
  - SVSCOREMEAN
- DGV: frequency in DGV, 8000 WGS

