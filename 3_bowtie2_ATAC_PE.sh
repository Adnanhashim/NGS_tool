#!/bin/bash

#SBATCH --account=nn8009k
# Project: test
#SBATCH --job-name=ATAC_PE
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --cpus-per-task=20
#SBATCH --output=/cluster/projects/nn8009k/Adnan/000_scripts/temp_slurm_output/3_bowtie2_PE_ATAC_%j.log
THREADS=20


set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

######### load modules #########
module --quiet purge
#module restore
module load deepTools/3.5.0-foss-2021a SAMtools/1.13-GCC-10.3.0  Bowtie2/2.4.4-GCC-10.3.0 BEDTools/2.30.0-GCC-10.3.0 BamTools/2.5.2-GCC-10.3.0  HTSlib/1.12-GCC-10.3.0 MultiQC/1.11-foss-2021a Subread/2.0.3-GCC-10.3.0 Trim_Galore/0.6.7-GCCcore-10.3.0 FastQC/0.11.9-Java-11

#module load Bowtie2/2.4.4-GCC-10.3.0
#module load SAMtools/1.12-GCC-10.3.0
#module load deepTools/3.3.1-intel-2019b-Python-3.7.4   


######### Inputs #########

fq1=$1
fq2=$2
#out=`basename $fq1 _R1_001_val_1.fq.gz`
out=`basename $fq1 _val_1.fq.gz`
DIR=/cluster/projects/nn8009k/Adnan
#READS=$3
OUTDIR=$3
specie=$4

###### indexes ############
#hg38_index="$DIR/00_Tools/Mads_indexes/hg38/bowtie2_canonical/hg38"
#hg19_index="$DIR/00_Tools/Mads_indexes/hg19/bowtie2_canonical/hg19"
#mm10_index="$DIR/00_Tools/Mads_indexes/mm10/bowtie2_canonical/mm10"

index=$DIR/00_Tools/Mads_indexes/${specie}/bowtie2_canonical/${specie}

####### effective sizes ######

if [ $specie = "mm10" ] ; then effective_genome_size="2652783500"; fi
if [ $specie = "hg38" ] ; then effective_genome_size="2913022398"; fi
if [ $specie = "hg19" ] ; then effective_genome_size="2864785220"; fi
if [ $specie = "mm9" ] ; then effective_genome_size="2620345972"; fi  
       

#echo "XXX"$mm10_index
#echo "YYY"$index
#echo "index is:" ${specie}_index
#echo "XXXXXXX:" $mm10_index
#echo "genome size is:" ${specie}_effective_genome_size

######### Outputs #########
 
#CHECK1="$READS/OUTPUT/3_bowtie2_PE_alignment/$out"
#CHECK1="$OUTDIR"

#if [ ! -e "$CHECK1" ] ; then
#    mkdir $OUTDIR
#    echo "output directory created :" $CHECK1 
#fi

ERR=$OUTDIR/$out.bowtie2_ATAC_PE.err
blacklisted_region=/cluster/projects/nn8009k/Adnan/00_Tools/blacklisted_regions/${specie}/*.bed

CHECK2="$OUTDIR/${out}_sorted.bam"
if [[ ! -f "$fq1" ]]; then echo "WRONG"; break; fi
if [[ ! -f "$fq2" ]]; then echo "WRONG"; break; fi


######### bowtie2 PE #########
date
    if [ ! -e "$CHECK2" ] ; then
        echo "Running bowtie2_PE, Sorting and indexing :" $out
        if bowtie2 -p $THREADS -x $index --very-sensitive -1 $fq1 -2 $fq2 1> $USERWORK/${out}.sam 2> $ERR
        #if bowtie2 -p $THREADS -x $mm10_index -1 $fq1 -2 $fq2 1> $USERWORK/${out}.sam 2> $ERR
           
           echo "Sorting"
           samtools view -bS $USERWORK/${out}.sam | samtools sort -@ $THREADS - -o $OUTDIR/${out}_tmp.bam

           echo "Removing blacklisted regions from alignment"
           bedtools intersect -v -abam $OUTDIR/${out}_tmp.bam -b $blacklisted_region | samtools sort -@ $THREADS - -o $OUTDIR/${out}_sorted.bam
           
           echo "indexing"
           samtools index -@ $THREADS $OUTDIR/${out}_sorted.bam $OUTDIR/${out}_sorted.bam.bai 
           
           ###### samtools Markdup
           echo "Running collate :" $out
           samtools collate -@ $THREADS -o $USERWORK/${out}_collate.bam $OUTDIR/${out}_sorted.bam

           echo "Running fixmate :" $out
           samtools fixmate -@ $THREADS -m $USERWORK/${out}_collate.bam $USERWORK/${out}_fixmate.bam

           echo "Running positionSort :" $out
           samtools sort -@ $THREADS -o $USERWORK/${out}_positionsort.bam $USERWORK/${out}_fixmate.bam

           echo "Running mardup :" $out
           samtools markdup -@ $THREADS -r -s $USERWORK/${out}_positionsort.bam $OUTDIR/${out}_dedup.bam
           samtools index $OUTDIR/${out}_dedup.bam $OUTDIR/${out}_dedup.bam.bai

           echo "creating bigwig"
           #bamCoverage -p $THREADS --bam $OUTDIR/${out}_dedup.bam -o $OUTDIR/${out}_dedup_RPGC.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize $effective_genome_size --ignoreForNormalization chrX chrM --ignoreDuplicates --extendReads
           bamCoverage -p $THREADS --bam $OUTDIR/${out}_dedup.bam -o $OUTDIR/${out}_dedup_RPKM.bw --binSize 10 --normalizeUsing RPKM --effectiveGenomeSize $effective_genome_size --ignoreForNormalization chrX chrM --ignoreDuplicates --extendReads
  
           rm $USERWORK/${out}.sam $OUTDIR/${out}_tmp.bam $USERWORK/${out}_collate.bam $USERWORK/${out}_fixmate.bam $USERWORK/${out}_positionsort.bam ;then

        echo Success
        else
            echo Failure
            exit
        fi    

        date
    fi    

echo "Done:" $out
    
