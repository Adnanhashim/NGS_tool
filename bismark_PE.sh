#!/bin/bash

#SBATCH --account=nn8009k
# Project: test
#SBATCH --job-name=bismark_PE
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --cpus-per-task=20
#SBATCH --output=/cluster/projects/nn8009k/Adnan/000_scripts/temp_slurm_output/bismark_%j.log
THREADS=20


set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

######### load modules #########
module --quiet purge
module restore
#module load Bowtie2/2.4.4-GCC-10.3.0
# module load SAMtools/1.12-GCC-10.3.0
#module load deepTools/3.3.1-intel-2019b-Python-3.7.4
#module --ignore-cache load deepTools/3.3.1-foss-2020a-Python-3.8.2  SAMtools/1.10-GCC-9.3.0 Bowtie2/2.4.1-GCC-9.3.0
#module load deepTools/3.3.1-foss-2020a-Python-3.8.2  SAMtools/1.10-GCC-9.3.0 Bowtie2/2.4.1-GCC-9.3.0
#module load HISAT2/2.1.0-foss-2018b
#module --ignore-cache load Bismark/0.23.1-GCC-10.3.0
module load Bismark/0.23.1-GCC-10.3.0 SAMtools/1.12-GCC-10.3.0
######### Inputs #########

fq1=$1
fq2=$2
#out=`basename $fq1 _R1_001_val_1.fq.gz`
out=`basename $fq1 _val_1.fq.gz`
DIR=/cluster/projects/nn8009k/Adnan
OUTDIR=$3
specie=$4

index=$DIR/00_Tools/bismark/${specie}

ERR=$OUTDIR/$out.bismark_PE.err

    CHECK2="$OUTDIR/${out}_sorted.bam"
        if [[ ! -f "$fq1" ]]; then echo "WRONG"; break; fi
        if [[ ! -f "$fq2" ]]; then echo "WRONG"; break; fi


######### bismark PE #########
date
    if [ ! -e "$CHECK2" ] ; then
        echo "-----Bismark-----Sample Name:"$out
            #if bismark --multicore $THREADS --bowtie2 --genome_folder $index -1 $fq1 -2 $fq2 --output_dir $OUTDIR/$out/ --temp_dir $USERWORK
            if bismark --multicore $THREADS --bowtie2 --genome_folder $index –non_directional –score_min L,0,-0.6, –un -1 $fq1 -2 $fq2 --output_dir $OUTDIR/$out/ --temp_dir $USERWORK #temporary 
                
                ##Deduplication
                    echo "-----deduplication-----Sample Name:"$out
                    echo "XXXXXXXXXXX input sample name:" $OUTDIR/$out/*.bam
                    echo "XXXXXXXXXXX output dir name:" $OUTDIR/$out/
                    echo "XXXXXXXXXXX output file name:" $out
                        deduplicate_bismark -p --bam --output_dir $OUTDIR/$out/ --outfile $out $OUTDIR/$out/*.bam

                ## methylation Extraction        
                    echo "-----methylation_extraction-----Sample Name:"$out
                        #bismark_methylation_extractor --multicore $THREADS --CX -p --ucsc -o $OUTDIR/$out/ --bedGraph $OUTDIR/$out/*deduplicated.bam
                        bismark_methylation_extractor --multicore $THREADS -p --ucsc -o $OUTDIR/$out/ --bedGraph $OUTDIR/$out/*deduplicated.bam #CpG
                        gunzip $OUTDIR/$out/*.cov.gz
                    
                  ##  echo "-----viewBS input-----Sample Name:"$out
                        #coverage2cytosine -CX -o $OUTDIR/$out/${out}_viewBS_input.cov --genome_folder $index $OUTDIR/$out/*.cov
                  ##      coverage2cytosine -o $OUTDIR/$out/${out}_viewBS_input.cov --genome_folder $index $OUTDIR/$out/*.cov #only CpGs
                
                ## Nucleotide coverage
                    echo "-----Bismark Nucleotide Coverage report-----Sample Name:"${out}
                        bam2nuc --genome_folder $index $OUTDIR/$out/*deduplicated.bam

                   # echo "-----Sorting by chromosome-----Sample Name:" ${out}
                   #     samtools sort -@ $THREADS -n $OUTDIR/$out/*deduplicated.bam -o $OUTDIR/$out/${out}_deduplicated_n_sorted.bam
           
                    echo "-----Bismark report-----Sample Name:"${out}
                        bismark2report ;then
                    echo Success
                else
                    echo Failure
                exit
            fi

        date
    fi    

echo "Done:" $out
