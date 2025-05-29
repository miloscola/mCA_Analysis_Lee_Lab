#!/bin/bash

############################################################################################################
#
#   Objective: Download test data and run the following scripts:
#		  Optional_step0_add_prefix_suffix.sh
#         step1_mCA_calling.sh      
#         optional_step2_mCA_compile.sh
#         step3_mCA_filter.R
#         Step4_allele_shift.sh
#
#   NOTE: line endings were changed from windows (CRLF) to UNIX (LF) for this file and all TopMed_WGS_mCA files.
#		  These files may not be properly altered on the github repo.
#
#   NOTE: this is not the optimal refrence genome for usage of this tool outside testing. Replace the file in step 1
#
#   USE: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/$gz_file
#
############################################################################################################
#   
#   Revision History: 
#
#   Date              Author               Title					                  Modification
#   05-25-2025        Mr.Milo Scola        Bioinformatics Programmer 3 Candidate      v1
#   
#   Contact								   Institute
#   miloscola@berkeley.edu				   UCSF
#
############################################################################################################
#
#   Directory Hierarchy
#
#   your_proj_working_directory
#   .
#   ├── down-sample                       ----> output BCFs after marker count-based downsampling
#   ├── GRCh38                            ----> prepared GRCh38 reference genome
#   ├── logs                              ----> to save SLURM output and error logs
#   │   └── *.out / *.err                 ----> logs for each SLURM job script
#   ├── mCA                               ----> mCA calling outputs and allele shift results
#   │   ├── AS                            ----> output folder for allele shift (AS) analysis results (e.g. Result.AS.ATM.bcf)
#   │   ├── gc                            ----> GC-annotated BCF files from mochatools
#   │   └── mCA_results                   ----> main folder for mCA calling results
#   │       ├── all.filtered.LoY.mCA.tsv  ----> final filtered Loss of Y (LoY) mCA table (step2_LoY output)
#   │       ├── all.mCA.tsv               ----> compiled mCA phenotypes across all samples (optional_step2 output)
#   │       ├── mCA_calls                 ----> individual mCA phenotype .tsv outputs from mocha
#   │       ├── mCA_stats                 ----> statistical summary tables output by mocha
#   │       └── mCA_vcf                   ----> mocha-generated VCF/BCF files with mCA calls
#   ├── mis                               ----> miscellaneous files, metadata, and sample lists
#   │   ├── GenomeWideHetCounts/          ----> raw heterozygosity count tables (before formatting)
#   │   ├── TopMed_sample_list            ----> plain list of sample IDs (used in prefix/suffix script)
#   │   ├── head.txt                      ----> template header file used for VCF reconstitution
#   │   ├── raw_VCF_list                  ----> output of optional_step0_add_prefix_suffix.sh
#   │   ├── vcf.list.txt                  ----> list of input VCF files for MoChA calling
#   │   └── *.het.count.*.list            ----> formatted sample+marker lists per race/gender (step1_select_sample)
#   ├── raw_data                          ----> input VCF/BCF raw data and modified headers
#   │   ├── cram_input                    ----> input .cra and .crai files for pre processing 
#   │   ├── *.bam                         ----> processed BAM files from CRAM 
#   │   ├── *.bcf                         ----> processed BCF files filtered and indexed
#   │   ├── head.*.txt                    ----> sample-specific VCF header files with SAMPLE replaced
#   │   └── head.*.vcf                    ----> concatenated header + body VCF files
#   └── scripts                           ----> bash and R scripts used throughout pipeline
#
#
#   Input: None
#
############################################################################################################

#TODO: ensure files are in correct place, may change WD
#TODO: add multiple files in testing phase

# Specify wd
host_data_dir="/data"
mkdir -p workspace
work_dir="/data"

# Make base directories
mkdir -p $work_dir/{logs,raw_data,mis,scripts,GRCh38,mCA}

## Download refrence genome if it does not already exist ##

# Set file names and directory
genome_ref_dir="$work_dir/GRCh38"
ref_file="GRCh38_full_analysis_set_plus_decoy_hla.fa"
gz_file="$ref_file.gz"


# Check if the unzipped .fna file exists
if [ -f "$genome_ref_dir/$ref_file" ]; then
  echo "Reference genome already exists at $genome_ref_dir/$ref_file. Skipping download."
else
  echo "Downloading and extracting GRCh38 reference genome..."
  wget -O "$genome_ref_dir/$ref_file" \
    ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
  gzip -d "$genome_ref_dir/$gz_file"
fi

# Generate index it doesn't exist
if [ ! -f "$genome_ref_dir/$ref_file.fai" ]; then
  echo "Reference index not found. Generating with samtools faidx..."
  samtools faidx "$ref_path"
else
  echo "Reference index already exists at ${ref_path}.fai. Skipping generation."
fi

## Download CRAM files if they do not already exist ##

# .cram/.crai input directory
cram_input_dir="$work_dir/raw_data/cram_input"
mkdir -p $cram_input_dir

# List of CRAM and CRAI file URLs
urls=(  
  "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/GBR/HG00114/alignment/HG00114.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage.cram"
  "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/GBR/HG00114/alignment/HG00114.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage.cram.crai"
)

# Loop through URLs
for url in "${urls[@]}"; do
  # Extract full filename from URL
  full_filename=$(basename "$url")

  # Extract just the sample ID from the filename (e.g., HG00096 from HG00096*.cram)
  sample_id=$(echo "$full_filename" | grep -oE 'HG[0-9]+')

  # Determine extension
  extension="${full_filename##*.}"

  # Set final filename
  filename="${sample_id}.${extension}"
  filepath="$cram_input_dir/$filename"

  #create temporary file
  tempfile="${filepath}.tmp"

  #if file is alrady there, do not download
  if [ -f "$filepath" ]; then
    echo "$filename already exists. Skipping download."
  else
    echo "Downloading $filename..."
    #if file is not downloaded properly, remove partial download 
    if wget -O "$tempfile" "$url"; then
      mv "$tempfile" "$filepath"
      echo "Downloaded and saved as $filename."
    else
      echo "Failed to download $filename. Removing partial file."
      rm -f "$tempfile"
    fi
  fi
done

echo "All CRAM and CRAI files are downloaded in $cram_input_dir."

## Create TopMed sample list ##

sample_list_ref_dir="$work_dir/mis"
echo -e "HG00114" > "$sample_list_ref_dir/TopMed_sample_list"

## Create VCF list ##

printf "HG00114.vcf" > "$sample_list_ref_dir/vcf.list.txt"

## generate VCF files ##

bash data/scripts/TopMed_WGS_mCA/optional_step-1_pre-process_CRAM.sh

## run optional_step0_add_prefix_suffix.sh ##

bash data/scripts/TopMed_WGS_mCA/optional_step0_add_prefix_suffix.sh

## run step1_mCA_calling.sh ##

bash data/scripts/TopMed_WGS_mCA/step1_mCA_calling.sh

## run optional_step2_mCA_compile.sh ##

bash data/scripts/TopMed_WGS_mCA/optional_step2_mCA_compile.sh

## run step3_mCA_filter ##

Rscript data/scripts/TopMed_WGS_mCA/step3_mCA_filter.R

## step4_allele_shift.sh ##
# NOTE: this step is not strictly necesissary and only 

#bash data/scripts/TopMed_WGS_mCA/step4_allele_shift.sh