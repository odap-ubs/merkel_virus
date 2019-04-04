#!/bin/bash
#$ -cwd
#$ -o /mnt/hydra/ubs/shared/users/Sandra/merkel_polyomavirus/res/log/
#$ -j y
#$ -S /bin/bash
# $ -pe smp 10
#$ -q all.q
export PATH=/share/apps/Perl/bin:$PATH
export LD_LIBRARY_PATH=../lib:${LD_LIBRARY_PATH}


### Folders
FASTQ="/mnt/hydra/ubs/shared/users/Sandra/IDIVAL"
TRIMMED="/mnt/hydra/ubs/shared/users/Sandra/merkel_polyomavirus/res/trimmed"
INDEX='/mnt/typhon/references/human/hg19/bowtie/hg19'
ALIGNMENTS='/mnt/hydra/ubs/shared/users/Sandra/merkel_polyomavirus/res/alignments'
UNMAPPED='/mnt/hydra/ubs/shared/users/Sandra/merkel_polyomavirus/res/alignments/unmapped'
VIRUS_ALIGN='/mnt/hydra/ubs/shared/users/Sandra/merkel_polyomavirus/res/virus_align'
VIRUS_MAPPED='/mnt/hydra/ubs/shared/users/Sandra/merkel_polyomavirus/res/virus_align/virus_mapped'
RES='/mnt/hydra/ubs/shared/users/Sandra/merkel_polyomavirus/res'

### Softwares
TRIMG="/share/apps/trim_galore"
CUTADAPT="/share/apps/cutadapt"
BOWTIE2='/share/apps/bowtie2'
SAMTOOLS='/share/apps/samtools/samtools'
BWA='/share/apps/bwa/bwa'


###########~~~~~~~~~~~~~~~~~~~~~~~Merkel polyomavirus project~~~~~~~~~~~~~~~~~~~~~#############
#################~~~~~~~~~~~~~~~~~~Finding virus in WES data~~~~~~~~~~~~~~~~~##################

### To look at colon tumor data, for seeing if there is possibly virus infection related to the 
### tumor incidence. We have WES data from colon but do not know if it works for finding virus
### integrated in the genome. So, before starting with colon data we do WES against the virus 
### from 15 tumor Merkel cell cancers WES data (4 of them have virus integrated in the genome)

### Sandra Garcia Mulero
### 21/08/2017



# --------------------------
# DOWNLOAD WES DATA (30 SAMPLES)
# -------------------------- 

cd /mnt/hydra/ubs/shared/users/Sandra/merkel_polyomavirus
ls -lhrt /mnt/hydra/ubs/shared/users/Sandra/IDIVAL | wc -l

# Info from the samples
cat /mnt/hydra/ubs/shared/users/Sandra/merkel_polyomavirus/MUESTRAS_MCC_SECUENCIADAS.csv
cat /mnt/hydra/ubs/shared/users/Sandra/merkel_polyomavirus/NGSCANT_08.csv

	# NUMBER OF SAMPLES IS 44

# --------------------------
# 1 FASTQC AND TRIMMING
# -------------------------- 

# 1.1 QC analysis

/share/apps/FastQC/fastqc /mnt/hydra/ubs/shared/users/Sandra/IDIVAL/C1728ACXX_2_1ss_1.fastq.gz
/share/apps/FastQC/fastqc /mnt/hydra/ubs/shared/users/Sandra/IDIVAL/C1VTNACXX_4_9ss_1.fastq.gz

firefox C1728ACXX_2_1ss_1_fastqc.html
firefox C1VTNACXX_4_9ss_1_fastqc.html

	# BAD QUALITY OF SEQUENCE READS QUALITY --> WE DO TRIMMING!!

# 1.2 Trimming the fastq


# tests:

# /share/apps/trim_galore/trim_galore --path_to_cutadapt /share/apps/cutadapt --quality 30 --paired --trim1 -o /mnt/hydra/ubs/shared/users/Sandra/merkel_polyomavirus/res/trimmed /mnt/hydra/ubs/shared/users/Sandra/IDIVAL/C1728ACXX_6_5ss_1.fastq.gz /mnt/hydra/ubs/shared/users/Sandra/IDIVAL/C1728ACXX_6_5ss_2.fastq.gz
# $TRIMG/trim_galore --path_to_cutadapt /share/apps/cutadapt --quality 30 --paired --trim1 -o $TRIMMED $FASTQ/C1728ACXX_6_5ss_1.fastq.gz $FASTQ/C1728ACXX_6_5ss_2.fastq.gz

# test worked!

for file in /mnt/hydra/ubs/shared/users/Sandra/IDIVAL/*_1.fastq.gz
do

base=$(basename $file "_1.fastq.gz")

$TRIMG/trim_galore --quality 30 --path_to_cutadapt $CUTADAPT --paired --trim1 -o $TRIMMED $FASTQ"/"$base"_1.fastq.gz" $FASTQ"/"$base"_2.fastq.gz"
echo trimming $base"_1.fastq.gz" $base"_2.fastq.gz" 

done

# --------------------------
# 2 ALIGNMENT: RUN BOWTIE2
# -------------------------- 

# 2.1 Generate genome indexes files

# BOWTIE2/bowtie2-build /mnt/typhon/references/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa REF/hg19.example

# El índice ya está creado y son estos files:
# hg19.1.bt2      hg19.2.bt2      hg19.3.bt2      hg19.4.bt2      hg19.rev.1.bt2  hg19.rev.2.bt2


# 2.2 Mapping reads to the genome

for file in /mnt/hydra/ubs/shared/users/Sandra/merkel_polyomavirus/res/trimmed/*_1_val_1.fq.gz
do

base=$(basename $file "_1_val_1.fq.gz")
echo aligning $base"_1_val_1.fq.gz" $base"_2_val_2.fq.gz"
$BOWTIE2/bowtie2 -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 -p 10 -x $INDEX -1 $TRIMMED/$base"_1_val_1.fq.gz" -2 $TRIMMED/$base"_2_val_2.fq.gz" | $SAMTOOLS view -bS - | $SAMTOOLS sort -m 1000000000 -O BAM -o $ALIGNMENTS/$base".sorted.bam"

done


# --------------------------
# 3 UNMAPPED READS 
# -------------------------- 

for file in /mnt/hydra/ubs/shared/users/Sandra/merkel_polyomavirus/res/alignments/*.bam
do

base=$(basename $file ".sorted.bam")
echo unmapped reads $base
$SAMTOOLS view -u -h -f 4 -F 264 $file > $UNMAPPED/tmps1.bam
$SAMTOOLS view -u -h -f 8 -F 260 $file > $UNMAPPED/tmps2.bam
$SAMTOOLS view -u -h -f 12 -F 256 $file > $UNMAPPED/tmps3.bam
$SAMTOOLS merge -u - tmps[123].bam | $SAMTOOLS sort -n - $UNMAPPED/$base"_unmapped_reads"
done

# --------------------------
# 4 ALIGN TO VIRAL GENOME
# -------------------------- 

for file in /mnt/hydra/ubs/shared/users/Sandra/merkel_polyomavirus/res/alignments/unmapped/*_unmapped_reads.bam
do

base=$(basename $file "_unmapped_reads.bam")
$BWA aln -b $VIRUS/Genome_Polyomavirus_MCC350.fasta $file > $VIRUS_ALIGN/$base".virus.sai"
$BWA samse $VIRUS/Genome_Polyomavirus_MCC350.fasta $VIRUS_ALIGN/$base".virus.sai" $file | $SAMTOOLS view - -bS -o $VIRUS_ALIGN/$base".virus.all.bam"
$SAMTOOLS view -@ 10 -b -F 0x04 $VIRUS_ALIGN/$base".virus.all.bam" -o $VIRUS_MAPPED/$base".virus.mapped.bam"
$SAMTOOLS view -@ 10 -b -f 0x04 $VIRUS_ALIGN/$base".virus.all.bam" -o $VIRUS_ALIGN/$base".virus.unmapped.bam"

done