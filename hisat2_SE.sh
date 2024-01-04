#!/bin/bash

#SBATCH --account=nn8009k
# Project: test
#SBATCH --job-name=hisat2_SE
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --cpus-per-task=20
#SBATCH --output=/cluster/projects/nn8009k/Adnan/000_scripts/temp_slurm_output/Hisat2_SE%j.log
THREADS=20


set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

######### load modules #########
module --quiet purge
#module restore
#module load Bowtie2/2.4.4-GCC-10.3.0
#module load SAMtools/1.12-GCC-10.3.0
#module load deepTools/3.3.1-intel-2019b-Python-3.7.4
#module --ignore-cache load deepTools/3.3.1-foss-2020a-Python-3.8.2  SAMtools/1.10-GCC-9.3.0 Bowtie2/2.4.1-GCC-9.3.0
module load HISAT2/2.2.1-gompi-2022a SAMtools/1.16.1-GCC-11.3.0


######### Inputs #########

fq1=$1
#fq2=$2
#out=`basename $fq1 _R1_001_val_1.fq.gz`
out=`basename $fq1 _trimmed.fq.gz`   
DIR=/cluster/projects/nn8009k/Adnan
OUTDIR=$2
specie=$3


#READS=$DIR/01_Madeleine/211027_A01447.A.Project_Fosslie-Libs18-2021-10-13/OUTPUT/2_trimglore/
#MAIN_OUT=$DIR/01_Madeleine/211027_A01447.A.Project_Fosslie-Libs18-2021-10-13/OUTPUT
#mm10_index=$DIR/00_Tools/Mads_indexes/mm10/bowtie2_canonical/mm10

index=$DIR/00_Tools/Hisat2_indexes/${specie}/genome

ERR=$OUTDIR/$out.hisat2_PE.err



CHECK2="$OUTDIR/${out}_sorted.bam"
#if [[ ! -f "$fq1" ]]; then echo "WRONG"; break; fi
#if [[ ! -f "$fq2" ]]; then echo "WRONG"; break; fi


######### Hisat2 SE #########
date
    if [ ! -e "$CHECK2" ] ; then
        echo "Running hisat2_SE, Sorting and indexing :" $out
        #if bowtie2 -p $THREADS -x $mm10_index -1 $fq1 -2 $fq2 1> ${out}.sam 2> $ERR  
        if hisat2 -p $THREADS -x $index --very-sensitive -U $fq1 -S $USERWORK/${out}.sam 2> $ERR
           
           echo "Sorting"
           samtools view -bS $USERWORK/${out}.sam | samtools sort -@ $THREADS - -o $OUTDIR/${out}_sorted.bam
           
           echo "indexing"
           samtools index -@ $THREADS $OUTDIR/${out}_sorted.bam $OUTDIR/${out}_sorted.bam.bai 
           
           rm $USERWORK/${out}.sam ;then
        echo Success
        else
            echo Failure
            exit
        fi    

        date
    fi    

echo "Done:" $out
