# NGS_tool
I have created NGS tool to analysis NGS data analyzing RNA-, ChIP-, ATAC- and bisulfite sequencing data on High-Performance Computing (HPC) and large-scale data storage to researchers in Norway Sigma2. https://www.sigma2.no/

## Introduction:


## Usage:

**Usage:**   sh NGS_tool.sh -i input_dir -o output_dir -m module -l layout -s Specie <br />
<br />
**Example:** sh NGS_tool.sh -i ./RawfastqFiles -o ./OUTPUT -m Alignment -l PE -s hg38  <br />
<br />
**Arguments:** <br />
 
**Parameter**           | **Description**
------------------------|-------------------------------------------------------------------
-i Input Directory      | A path of input directory
-o Output Directory     | A path of output directory
-m module               | module you want to run, e.g., fastqc, Trimglore, Alignment (ChIP), Hisat2 (RNAseq), bismark (WGBS), ATAC"
-l layout               | layout of reads, e.g., SE, PE
-s Specie               | select the species, e.g., hg38, hg19, mm10, mm9, bosTau8, susScr11
-h/--help               | Show this help message and exit

## fastq file name:

fastq file name should end with _R1_001.fastq.gz and _R2_001.fastq.gz you can following command to rename your fastq files

**PE:**

for i in *_1.fastq.gz; do   mv "$i" "`echo $i | sed "s/_1.fastq.gz/_R1_001.fastq.gz/"`"; done
for i in *_2.fastq.gz; do   mv "$i" "`echo $i | sed "s/_2.fastq.gz/_R2_001.fastq.gz/"`"; done

for i in *_1.fq.gz; do mv "$i" "$(echo $i | sed 's/_1.fq.gz/_R1_001.fastq.gz/')"; done
for i in *_2.fq.gz; do mv "$i" "$(echo $i | sed 's/_2.fq.gz/_R2_001.fastq.gz/')"; done



**SE:**

for i in *gz; do   mv "$i" "`echo $i | sed "s/\.fastq.gz/_R1_001.fastq.gz/"`"; done
