#!/bin/bash

############################################################################################################
#
#   Objective: pre-process CRAM files to get VCFs         
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
#   ├── GRCh38              ----> to save the prepared GRCh38 resources 
#   ├── scripts 			----> to save the bash scripts  
#   ├── mCA             	----> to save the mChA calling result files
#   |── raw_data            ----> to store the input vcf/bcf raw data 
#   │   └── cram_input      ----> input .cra and .crai files for pre processing 
#	|── logs				----> to save logs for checkup
#   ├── mis    				----> to save sample list, metadata and miscellaneous		
#
#
#   Input: A CRAM sample list .txt file, each row should be a single sample.
#   
#   Like below:
#   			NWD123654
#   			NWD456987
#   			NWD159753
#   			......
#   			......
#   
#   A set of paired cram/crai files
#
############################################################################################################
echo -e "script:optional_step-1_pre-process_CRAM started at $(date)/n"

# specify your wd
wd="/data"
echo -e "Working Directory: $wd\n"

# Make sure the directories are ready
mkdir -p $wd/{logs,raw_data,mis,scripts,GRCh38,mCA}

# Make sure the TopMed sample list exists, it should be one row one sample ID
sample_dir="$wd/mis"

# Path to GRCh38 reference genome
ref="$wd/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

# Ensure reference genome index exists
if [ ! -f "$ref.fai" ]; then
    echo "Reference index not found, generating..."
    samtools faidx $ref
fi

# .cram/.crai input directory
cram_input_dir="$wd/raw_data/cram_input"
mkdir -p $cram_input_dir

# Loop over each sample ID in TopMed_sample_list
while read sample_id; do
    echo -e "Processing sample: $sample_id\n"

    # Define input CRAM path
    cram_file="$cram_input_dir/${sample_id}.cram"

    # Convert CRAM to BAM
    bam_file="$wd/raw_data/${sample_id}.bam"
    samtools view -T $ref -b -o $bam_file $cram_file

    # Sort BAM
    sorted_bam="$wd/raw_data/${sample_id}.sorted.bam"
    samtools sort -o $sorted_bam $bam_file
    rm $bam_file  # optional cleanup

    # Index sorted BAM
    samtools index $sorted_bam

    # Generate VCF with GT, AD, DP (MoChA expects these fields)
    vcf_out="$wd/raw_data/${sample_id}.vcf"
    bcftools mpileup -Ou -f $ref -a FORMAT/AD,FORMAT/DP $sorted_bam | \
    bcftools call -m -Ov -o $vcf_out

    echo -e "Sample $sample_id processed to VCF: $vcf_out\n"

done < "$sample_dir/TopMed_sample_list"

echo -e "script:optional_step-1_pre-process_CRAM ended at $(date)\n"