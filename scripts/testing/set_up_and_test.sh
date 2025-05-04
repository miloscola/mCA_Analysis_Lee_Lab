
############################################################################################################
#
#   Objective: Download test data and run the following scripts:
#		  Optional_step0_add_prefix_suffix.sh
#         step1_mCA_calling.sh      
#         optional_step2_mCA_compile.sh
#         step3_mCA_filter.R
#         Step4_allele_shift.sh
#
#   NOTE: line endings were changed from windows (CRLF) to UNIX (LF) for this file and all TopMed_WGS_mCA files
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
#   ├── GRCh38                            ----> prepared GRCh38 reference genome
#   ├── scripts                           ----> bash and R scripts used throughout pipeline
#   ├── mCA                               ----> mCA calling outputs and allele shift results
#   │   ├── gc                            ----> GC-annotated BCF files from mochatools
#   │   ├── mCA_results                   ----> main folder for mCA calling results
#   │   │   ├── mCA_calls                 ----> individual mCA phenotype .tsv outputs from mocha
#   │   │   ├── mCA_stats                 ----> statistical summary tables output by mocha
#   │   │   ├── mCA_vcf                   ----> mocha-generated VCF/BCF files with mCA calls
#   │   │   ├── all.mCA.tsv               ----> compiled mCA phenotypes across all samples (optional_step2 output)
#   │   │   └── all.filtered.LoY.mCA.tsv  ----> final filtered Loss of Y (LoY) mCA table (step2_LoY output)
#   │   └── AS                            ----> output folder for allele shift (AS) analysis results (e.g. Result.AS.ATM.bcf)
#   ├── raw_data                          ----> input VCF/BCF raw data and modified headers
#   │   ├── head.*.txt                    ----> sample-specific VCF header files with SAMPLE replaced
#   │   ├── head.*.vcf                    ----> concatenated header + body VCF files
#   │   └── *.bcf                         ----> processed BCF files filtered and indexed
#   ├── down-sample                       ----> output BCFs after marker count-based downsampling
#   ├── mis                               ----> miscellaneous files, metadata, and sample lists
#   │   ├── vcf.list.txt                  ----> list of input VCF files for MoChA calling
#   │   ├── TopMed_sample_list            ----> plain list of sample IDs (used in prefix/suffix script)
#   │   ├── raw_VCF_list                  ----> output of optional_step0_add_prefix_suffix.sh
#   │   ├── head.txt                      ----> template header file used for VCF reconstitution
#   │   ├── *.het.count.*.list            ----> formatted sample+marker lists per race/gender (step1_select_sample)
#   │   └── GenomeWideHetCounts/          ----> raw heterozygosity count tables (before formatting)
#   ├── logs                              ----> to save SLURM output and error logs
#   │   └── *.out / *.err                 ----> logs for each SLURM job script
#
#
#   Input: None
#
############################################################################################################

# Specify wd
wd=$(basename "$(dirname "$(dirname "$PWD")")")
echo -e "Working Directory: $wd\n"

# Make base directories
mkdir -p $wd/{logs,raw_data,mis,scripts,GRCh38,mCA}

# Download refrence genome
wget -O- ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz 
gzip -d $wd/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
samtools faidx $wd/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

