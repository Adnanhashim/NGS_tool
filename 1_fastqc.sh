#!/bin/bash

#SBATCH --account=nn8009k
#SBATCH --job-name=fastqc
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/cluster/projects/nn8009k/Adnan/000_scripts/temp_slurm_output/1_fastqc_%j.log

set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

# Load required modules
module purge
module load deepTools/3.5.0-foss-2021a SAMtools/1.13-GCC-10.3.0 Bowtie2/2.4.4-GCC-10.3.0 BEDTools/2.30.0-GCC-10.3.0 BamTools/2.5.2-GCC-10.3.0 HTSlib/1.12-GCC-10.3.0 MultiQC/1.11-foss-2021a Subread/2.0.3-GCC-10.3.0 Trim_Galore/0.6.7-GCCcore-10.3.0 FastQC/0.11.9-Java-11

# Input file and base name
FILE="$1"
BASE=$(basename "$FILE" .fastq.gz)

# Output directory
OUTDIR="$2"

# Define the path to the FastQC executable
FASTQC=$(which fastqc)

# Create an error log file
FASTQC_ERR="$OUTDIR/$BASE.fastqc.err"

# Print a message to the log
echo "Running FastQC: $BASE"

# Run FastQC
$FASTQC --outdir "$OUTDIR" "$FILE" 2> "$FASTQC_ERR"

# Check if FastQC completed successfully
if [ $? -eq 0 ]; then
    echo "FastQC completed successfully for: $BASE"
else
    echo "FastQC encountered errors. Check the log: $FASTQC_ERR"
fi


