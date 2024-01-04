#!/bin/bash

#SBATCH --account=nn8009k
# Project: test
#SBATCH --job-name=bowtie2_SE
#SBATCH --time=16:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --cpus-per-task=20
#SBATCH --output=/cluster/projects/nn8009k/Adnan/000_scripts/temp_slurm_output/3_bowtie2_SE_%j.log
THREADS=20


set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

######### load modules #########
module --quiet purge
module load deepTools/3.5.0-foss-2021a SAMtools/1.13-GCC-10.3.0  Bowtie2/2.4.4-GCC-10.3.0 BEDTools/2.30.0-GCC-10.3.0 BamTools/2.5.2-GCC-10.3.0  HTSlib/1.12-GCC-10.3.0 MultiQC/1.11-foss-2021a Subread/2.0.3-GCC-10.3.0 Trim_Galore/0.6.7-GCCcore-10.3.0 FastQC/0.11.9-Java-11


######### Inputs #########

fq1=$1
#fq2=$2
#out=`basename $fq1 _R1_001_val_1.fq.gz`
out=`basename $fq1 _trimmed.fq.gz`
DIR=/cluster/projects/nn8009k/Adnan
#READS=$2
OUTDIR=$2
specie=$3

###### indexes ############
#index=$DIR/00_Tools/Mads_indexes/${specie}/bowtie2_canonical/${specie}
index=$DIR/00_Tools/Mads_indexes/${specie}/bowtie2_canonical/${specie}
#index=$DIR/00_Tools/Mads_indexes/mm10/bowtie2_canonical/mm10

####### effective sizes ######

if [ $specie = "mm10" ] ; then effective_genome_size="2652783500"; fi
if [ $specie = "hg38" ] ; then effective_genome_size="2913022398"; fi
if [ $specie = "hg19" ] ; then effective_genome_size="2864785220"; fi
if [ $specie = "mm9" ] ; then effective_genome_size="2620345972"; fi  
       
######### Outputs #########
 

ERR=$OUTDIR/$out.bowtie2_SE.err
#blacklisted_region=/cluster/projects/nn8009k/Adnan/00_Tools/blacklisted_regions/${specie}/*.bed
blacklisted_region=/cluster/projects/nn8009k/Adnan/00_Tools/blacklisted_regions/mm10/*.bed

CHECK2="$OUTDIR/${out}_sorted.bam"

######### bowtie2 SE #########
date
    if [ ! -e "$CHECK2" ] ; then
        echo "Running bowtie2_SE, Sorting and indexing :" $out
        
        if bowtie2 -p $THREADS -x $index -U $fq1 -S $USERWORK/${out}.sam 2> $ERR
           
           echo "Sorting"
           samtools view -bhS $USERWORK/${out}.sam | samtools sort -@ $THREADS - -o $OUTDIR/${out}_tmp.bam

           echo "Removing blacklisted regions from alignment"
           bedtools intersect -v -abam $OUTDIR/${out}_tmp.bam -b $blacklisted_region | samtools sort -@ $THREADS - -o $OUTDIR/${out}_sorted.bam
           
           echo "indexing"
           samtools index -@ $THREADS $OUTDIR/${out}_sorted.bam $OUTDIR/${out}_sorted.bam.bai 

           #echo "creating bigwig"
           #bamCoverage -p $THREADS --bam $OUTDIR/${out}_sorted.bam -o $OUTDIR/${out}_RPKM.bw --binSize 10 --normalizeUsing RPKM --effectiveGenomeSize $effective_genome_size --ignoreForNormalization chrX chrM --ignoreDuplicates

           echo "Running mardup :" $out
           samtools markdup -@ $THREADS -r -s $OUTDIR/${out}_sorted.bam $OUTDIR/${out}_dedup.bam
           samtools index $OUTDIR/${out}_dedup.bam $OUTDIR/${out}_dedup.bam.bai
           
           echo "creating bigwig"
           bamCoverage -p $THREADS --bam $OUTDIR/${out}_dedup.bam -o $OUTDIR/${out}_dedup_RPKM.bw --binSize 10 --normalizeUsing RPKM --effectiveGenomeSize $effective_genome_size --ignoreForNormalization chrX chrM --ignoreDuplicates
           rm $USERWORK/${out}.sam $OUTDIR/${out}_tmp.bam;then
            
        echo Success
        else
            echo Failure
            exit
        fi    

        date
    fi    

echo "Done:" $out
