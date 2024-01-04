#!/bin/bash

#SBATCH --account=nn8009k
#SBATCH --job-name=trimglore
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=/cluster/projects/nn8009k/Adnan/000_scripts/temp_slurm_output/2_trimglore_%j.log

set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

# Load necessary modules
module --quiet purge
module load deepTools/3.5.0-foss-2021a SAMtools/1.13-GCC-10.3.0 Bowtie2/2.4.4-GCC-10.3.0 BEDTools/2.30.0-GCC-10.3.0 BamTools/2.5.2-GCC-10.3.0 HTSlib/1.12-GCC-10.3.0 MultiQC/1.11-foss-2021a Subread/2.0.3-GCC-10.3.0 Trim_Galore/0.6.7-GCCcore-10.3.0 FastQC/0.11.9-Java-11

# Define input variables
fq1="$1"
fq2="$2"
OUTDIR="$3"
out="$4"
DIR="/cluster/projects/nn8009k/Adnan"
ERR="$OUTDIR/$out.trimglore.err"

# Check if input files exist
if [[ ! -f "$fq1" ]] || [[ ! -f "$fq2" ]]; then
    echo "Error: Input files do not exist."
    exit 1
fi

# Check for output files
CHECK2="${OUTDIR}/${out}*_val_*.fq.gz"

# Only proceed if the output files don't exist
if [[ ! -e "$CHECK2" ]]; then
    echo "Running Trimglore: $out"
    if trim_galore -j 4 --nextseq=20 --paired "$fq1" "$fq2" -o "$OUTDIR" --fastqc_args "--outdir $OUTDIR/fastqc" 2> "$ERR"; then
        echo "Trimglore completed successfully."
    else
        echo "Trimglore failed."
        exit 1
    fi

    date
fi

echo "Done: $out"

