#!/bin/bash

# Define job parameters
#SBATCH --account=nn8009k
#SBATCH --job-name=bowtie2_PE
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --cpus-per-task=20
#SBATCH --output=/cluster/projects/nn8009k/Adnan/000_scripts/temp_slurm_output/3_bowtie2_PE_%j.log

# Set thread count
THREADS=20

# Enable error checking
set -o errexit
set -o nounset

# Load required modules
module --quiet purge
module restore
module load deepTools/3.5.0-foss-2021a SAMtools/1.13-GCC-10.3.0  Bowtie2/2.4.4-GCC-10.3.0 BEDTools/2.30.0-GCC-10.3.0 BamTools/2.5.2-GCC-10.3.0  HTSlib/1.12-GCC-10.3.0 MultiQC/1.11-foss-2021a Subread/2.0.3-GCC-10.3.0 Trim_Galore/0.6.7-GCCcore-10.3.0 FastQC/0.11.9-Java-11

# Define input variables
fq1=$1
fq2=$2
out=$(basename "$fq1" _val_1.fq.gz)
DIR="/cluster/projects/nn8009k/Adnan"
OUTDIR="$3"
specie="$4"

# Define indexes
index="$DIR/00_Tools/Mads_indexes/${specie}/bowtie2_canonical/${specie}"

# Define effective genome sizes
case $specie in
  "mm10") effective_genome_size="2652783500";;
  "hg38") effective_genome_size="2913022398";;
  "hg19") effective_genome_size="2864785220";;
  "mm9") effective_genome_size="2620345972";;
  "bosTau8") effective_genome_size="2640168518";;
  "susScr11") effective_genome_size="2435278676";;
  "TC1") effective_genome_size="2699493483";;
  *) echo "Unsupported species: $specie"; exit 1;;
esac

# Define blacklisted region file
#blacklisted_region="/cluster/projects/nn8009k/Adnan/00_Tools/blacklisted_regions/${specie}/*.bed"
#blacklisted_region=$(find /cluster/projects/nn8009k/Adnan/00_Tools/blacklisted_regions/${specie} -name '*.bed' -type f)
#blacklisted_region=/cluster/projects/nn8009k/Adnan/00_Tools/blacklisted_regions/${specie}/*.bed

# Define output paths
CHECK2="$OUTDIR/${out}_sorted.bam"
ERR="$OUTDIR/${out}.bowtie2_PE.err"

# Check if input files exist
if [[ ! -f "$fq1" || ! -f "$fq2" ]]; then
  echo "Input files do not exist."
  exit 1
fi

# Check if output directory exists; create it if not
if [ ! -e "$OUTDIR" ]; then
  mkdir -p "$OUTDIR"
  echo "Output directory created: $OUTDIR"
fi


# Run bowtie2_PE and processing steps
if [[ ! -e "$CHECK2" ]]; then
  echo "Running bowtie2_PE, Sorting and indexing: $out"
  bowtie2 -p "$THREADS" -x "$index" -1 "$fq1" -2 "$fq2" 1> "$USERWORK/${out}.sam" 2> "$ERR"
  
  echo "Sorting :" $out
  samtools view -bhS "$USERWORK/${out}.sam" | samtools sort -@ "$THREADS" - -o "$OUTDIR/${out}_tmp.bam"

  # Optionally remove blacklisted regions based on the species
  if [[ "$specie" != "susScr11" && "$specie" != "TC1" ]]; then
    blacklisted_region=$(find /cluster/projects/nn8009k/Adnan/00_Tools/blacklisted_regions/${specie} -name '*.bed' -type f)
    bedtools intersect -v -abam "$OUTDIR/${out}_tmp.bam" -b "$blacklisted_region" | samtools sort -@ "$THREADS" - -o "$OUTDIR/${out}_sorted.bam"
  else
    mv "$OUTDIR/${out}_tmp.bam" "$OUTDIR/${out}_sorted.bam"
  fi

  echo "Indexing :" $out
  samtools index -@ "$THREADS" "$OUTDIR/${out}_sorted.bam" "$OUTDIR/${out}_sorted.bam.bai"

  # Perform other processing steps (collate, fixmate, positionSort, markdup, create bigwig)

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
                bamCoverage -p $THREADS --bam $OUTDIR/${out}_dedup.bam -o $OUTDIR/${out}_dedup_RPKM.bw --binSize 10 --normalizeUsing RPKM --effectiveGenomeSize $effective_genome_size --ignoreForNormalization chrX chrM --ignoreDuplicates --extendReads


  rm "$USERWORK/${out}.sam" "$OUTDIR/${out}_tmp.bam" "$USERWORK/${out}_collate.bam" "$USERWORK/${out}_fixmate.bam" "$USERWORK/${out}_positionsort.bam"

  echo "Success"
else
  echo "Output file already exists: $CHECK2"
fi


echo "Done: $out"
