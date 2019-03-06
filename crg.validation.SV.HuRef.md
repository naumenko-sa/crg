# HuRef SVs

Benchmark (DEL calls): http://www.cell.com/action/showImagesData?pii=S0002-9297%2817%2930496-2
Raw data: https://www.ncbi.nlm.nih.gov/sra/SRX5395595[accn]

## Prepare input data:
```
module load sratoolkit
prefetch -c SRR8595488 --max-size 60GB
qsub ~/bioscripts/fastq.sra2fq.sh -v srr=SRR8595488,sample=huref_blood1
```
Result
```
38G	huref_blood1_1.fq.gz
42G	huref_blood1_2.fq.gz
```
