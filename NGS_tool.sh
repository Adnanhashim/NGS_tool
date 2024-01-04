#!/bin/sh

# Source the functions script to make the functions available
source "./functions.sh"

# Define the help function
helpFunction() {
  echo ""
  echo "Usage: $0 -i 'Input_dir' -o 'Output_dir' -m 'Module' -l 'layout' -s 'specie'"
  echo -e "\t-i Input Directory"
  echo -e "\t-o Output Directory"
  echo -e "\t-m module you want to run, e.g., fastqc, Trimglore, Alignment (ChIP), Hisat2 (RNAseq), bismark (WGBS), ATAC"
  echo -e "\t-l layout of reads, e.g., SE, PE"
  echo -e "\t-s select the species, e.g., hg38, hg19, mm10, mm9, bosTau8, susScr11, TC1"
  exit 1 # Exit script after printing help
}

# Parse command line options and arguments
while getopts "i:o:m:l:s:" opt; do
  case "$opt" in
    i ) Input_dir="$OPTARG" ;;
    o ) Output_dir="$OPTARG" ;;
    m ) Module="$OPTARG" ;;
    l ) layout="$OPTARG" ;;
    s ) specie="$OPTARG" ;;
    ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
  esac
done

# Check if mandatory parameters are empty and display help if necessary
if [ -z "$Input_dir" ] || [ -z "$Output_dir" ] || [ -z "$Module" ]; then
  echo "Some or all of the parameters are empty: Input_dir, Output_dir, and module are mandatory"
  helpFunction
fi

# Execute the selected analysis module
case "$Module" in
  "fastqc")
    run_fastqc
    ;;
  "Trimglore")
    run_trimglore
    ;;
  "Alignment")
    run_alignment
    ;;
  "Hisat2")
    run_hisat2
    ;;
  "bismark")
    run_bismark
    ;;
  "ATAC")
    run_atac
    ;;
  # Add cases for other modules as needed
  *)
    echo "Please select a valid module."
    exit 1
    ;;
esac









































