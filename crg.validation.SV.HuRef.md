# Validation of structural variant (SV) calling in bcbio using HuRef benchmark

Benchmark (DEL calls): http://www.cell.com/action/showImagesData?pii=S0002-9297%2817%2930496-2\ 
Raw data: https://www.ncbi.nlm.nih.gov/sra/SRX5395595[accn]

## 1. Prepare input data:
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

## 2. Align to grch37 with decoy

```
details:
- algorithm:
    aligner: bwa
    effects: false
    mark_duplicates: true
    realign: false
    recalibrate: false
    save_diskspace: true
    tools_on:
    - svplots
    - qualimap
    variantcaller: false
  analysis: variant2
  description: huref_blood1
  files:
  - /path/huref/input/huref_blood1_1.fq.gz
  - /path/huref/input/huref_blood1_2.fq.gz
  genome_build: GRCh37d5
  metadata:
    batch: huref
fc_name: huref
resources:
  default:
    cores: 7
    jvm_opts:
    - -Xms750m
    - -Xmx7000m
    memory: 7G
upload:
  dir: ../final
```

## 3. Filter out decoy reads

[cre.bam.remove_decoy_reads.sh](https://github.com/naumenko-sa/cre/blob/master/cre.bam.remove_decoy_reads.sh)

## 4. Run SV calling (all tools, or each tool in a separate project)

```

```