# NGS_tool
I have created NGS tool to analysis NGS data analyzing RNA-, ChIP-, ATAC- and bisulfite sequencing data on High-Performance Computing (HPC) and large-scale data storage to researchers in Norway Sigma2. https://www.sigma2.no/

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

