# Detection of Merkel cell Polyomavirus from WES data
A pipeline for detetction of Merkel Cell Polyomavirus insertions in human DNA by WES analysis

## Introduction
Merkel cell carcinoma (MCC) is highly malignant neuroendocrine tumor of the skin arising from the mechanoreceptor Merker cells. It is caused by the Merkel cell Polyomavirus (MCV), a DNA virus that is detected in 75–89% of MCC. Normally, MCV insertions in MCC genomes were detected by PCR. Here we shown a 2-alignment-steps methodology to identify MCV+ tumors from whole exome sequencing (WES) data and using the MCV genome.

## Requeriments
Written in Bash.
Requires:
* TrimGalore
* Cutadapt
* Bowtie2
* Samtools
* BWA

## Steps
### 1. Quality pre-process
```
# FastQC
fastqc C1VTNACXX_4_9ss_1.fastq.gz
fastqc C1VTNACXX_4_9ss_2.fastq.gz

firefox C1VTNACXX_4_9ss_1_fastqc.html
firefox C1VTNACXX_4_9ss_2.fastq.html

# Trimming to remove bad quality 
for file in /mnt/hydra/ubs/shared/users/Sandra/IDIVAL/*_1.fastq.gz; do 
  base=$(basename $file "_1.fastq.gz")
  $TRIMG/trim_galore --quality 30 --path_to_cutadapt $CUTADAPT --paired --trim1 -o $TRIMMED $FASTQ"/"$base"_1.fastq.gz" $FASTQ"/"$base"_2.fastq.gz"
  echo trimming $base"_1.fastq.gz" $base"_2.fastq.gz" 
done
```
### 2.  
### 3. 
