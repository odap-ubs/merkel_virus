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
The Raw data from the WES fastq files are first trimmed if necessary with TrimGalore (v0.4.0).
```bash
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
### 2.  Alignemnt 1 - against human
Then, trimmed reads are aligned to the reference human genome build hg19/GRCh37, previously indexed, using Bowtie 2.0 (v2.2.5) ─alignment 1─. 
```bash
# 2.1 Generate genome indexes files
$BOWTIE2/bowtie2-build /mnt/typhon/references/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa REF/hg19.example
# 2.2 Mapping reads to the genome
for file in /mnt/hydra/ubs/shared/users/Sandra/merkel_polyomavirus/res/trimmed/*_1_val_1.fq.gz; do
  base=$(basename $file "_1_val_1.fq.gz")
  echo aligning $base"_1_val_1.fq.gz" $base"_2_val_2.fq.gz"
  $BOWTIE2/bowtie2 -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 -p 10 -x $INDEX -1 $TRIMMED/$base"_1_val_1.fq.gz" -2 $TRIMMED/$base"_2_val_2.fq.gz" | $SAMTOOLS view -bS - | $SAMTOOLS sort -m 1000000000 -O BAM -o $ALIGNMENTS/$base".sorted.bam"
done
```
### 3. Retrieve unmapped reads
From this first alignment, unmapped reads are selected with Samtools (v1.3.1).
```bash
for file in /mnt/hydra/ubs/shared/users/Sandra/merkel_polyomavirus/res/alignments/*.bam; do
  base=$(basename $file ".sorted.bam")
  echo unmapped reads $base
  # 1- An unmapped read whose mate is mapped
  $SAMTOOLS view -u -h -f 4 -F 264 $file > $UNMAPPED/tmps1.bam
  # 2- A mapped read who’s mate is unmapped
  $SAMTOOLS view -u -h -f 8 -F 260 $file > $UNMAPPED/tmps2.bam
  # 3- Both reads of the pair are unmapped
  $SAMTOOLS view -u -h -f 12 -F 256 $file > $UNMAPPED/tmps3.bam
  # 4- Merge
  $SAMTOOLS merge -u - tmps[123].bam | $SAMTOOLS sort -n - $UNMAPPED/$base"_unmapped_reads"
done
```
### 4- Alignment 2 - against viral genome
```bash
for file in /mnt/hydra/ubs/shared/users/Sandra/merkel_polyomavirus/res/alignments/unmapped/*_unmapped_reads.bam; do
  base=$(basename $file "_unmapped_reads.bam")
  $BWA aln -b $VIRUS/Genome_Polyomavirus_MCC350.fasta $file > $VIRUS_ALIGN/$base".virus.sai"
  $BWA samse $VIRUS/Genome_Polyomavirus_MCC350.fasta $VIRUS_ALIGN/$base".virus.sai" $file | $SAMTOOLS view - -bS -o $VIRUS_ALIGN/$base".virus.all.bam"
  $SAMTOOLS view -@ 10 -b -F 0x04 $VIRUS_ALIGN/$base".virus.all.bam" -o $VIRUS_MAPPED/$base".virus.mapped.bam"
  $SAMTOOLS view -@ 10 -b -f 0x04 $VIRUS_ALIGN/$base".virus.all.bam" -o $VIRUS_ALIGN/$base".virus.unmapped.bam"
done
``
