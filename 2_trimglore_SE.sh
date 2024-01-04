#!/bin/bash

#SBATCH --account=nn8009k
#SBATCH --job-name=trimglore
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=/cluster/projects/nn8009k/Adnan/000_scripts/temp_slurm_output/2_trimglore_%j.log

#######SBATCH --output=//cluster/projects/nn8009k/Adnan/01_Madeleine/211027_A01447.A.Project_Fosslie-Libs18-2021-10-13/OUTPUT/0_slurm_outputs/2_trimglore_%j.log


set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

######### load modules #########
module --quiet purge
module load deepTools/3.5.0-foss-2021a SAMtools/1.13-GCC-10.3.0  Bowtie2/2.4.4-GCC-10.3.0 BEDTools/2.30.0-GCC-10.3.0 BamTools/2.5.2-GCC-10.3.0  HTSlib/1.12-GCC-10.3.0 MultiQC/1.11-foss-2021a Subread/2.0.3-GCC-10.3.0 Trim_Galore/0.6.7-GCCcore-10.3.0 FastQC/0.11.9-Java-11


######### Inputs #########

fq1=$1
out=`basename $fq1 .fastq.gz`
DIR=/cluster/projects/nn8009k/Adnan


######### Outputs ######### 
OUTDIR=$2

if [ ! -e "$OUTDIR" ] ; then
    mkdir -p "$OUTDIR"
fi

ERR=$OUTDIR/$out.trimglore.err

CHECK2="${out}*_trimmed.fq.gz"


######### Trim Galore #########
date
    if [ ! -e "$CHECK2" ] ; then
        echo "Running Trimglore :" $out
        if trim_galore -j 4 $fq1 -o $OUTDIR --fastqc_args "--outdir $OUTDIR/fastqc" 2> $ERR ;then  #no polyG issue
        #if trim_galore -j 4 --nextseq=20 --paired $fq1 $fq2 -o $OUTDIR --fastqc_args "--outdir $READS/OUTPUT/2_trimglore/fastqc" 2> $ERR ;then  #nextseq polyG issue
            echo Success
        else
            echo Failure
            exit
        fi    

        date
    fi    

echo "Done:" $out
