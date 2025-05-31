# TOPMed mCA Pipeline

This repository contains a pipeline for detecting mosaic chromosomal alterations (mCAs) in whole-genome sequencing data. Scripts for detecting loss of the Y chromosome (LoY) and allele shifts are included but not used in the test pipeline. The testing script limits its analysis to chromosome 22 for expediency.

This pipeline reproduces the analysis described in:

> **Terao, C., Wang, B., Hirata, J. et al. Mosaic chromosomal alterations across 756,000 genomes in the TOPMed program.** *Nature Genetics* **55**, 1903–1915 (2023). [https://doi.org/10.1038/s41588-023-01553-1](https://doi.org/10.1038/s41588-023-01553-1)

The pipeline consists of several modular steps:

### 1. **optional_step-1_pre-process_CRAM**
- **Function:** Generates properly formatted VCF files for MoCHA from CRAM inputs
- **Input:** `.cram` and `.crai` files, reference genome (`GRCh38`).
- **Output:** Preprocessed VCF files.

### 2. **optional_step0_add_prefix_suffix**
- **Function:** Standardizes sample IDs and filenames to ensure downstream compatibility.
- **Input:** Sample list and metadata (`mis/TopMed_sample_list`).
- **Output:** Normalized file names and sample map.

### 3. **step1_mCA_calling**
- **Function:** Runs variant calling, applies depth and GC filters, and performs mCA detection using [MOCHA](https://github.com/freeseek/mocha).
- **Input:** VCF files, VCF list .txt file, and reference genome.
- **Output:** Per-sample `.bcf` with mCA annotations, summary of results in mca/

### 4. **optional_step2_mCA_compile**
- **Function:** Aggregates per-sample mCA calls into a unified tabular file.
- **Output:** `proj1.all.mCA.tsv` summarizing all calls.

## Test Data

The provided test dataset includes a publicly available CRAM file:

- **Sample 1:** HG00234 from the 1000 Genomes Project.
- **Sample 2:** HG00114 from the 1000 Genomes Project.
- **Sample 3:** HG00096 from the 1000 Genomes Project.

**Expected Result:** These samples are a healthy controls; **no mCAs** are expected to be detected. Public data with mCAs available without licence could not be found.

This dataset is used to validate that the pipeline runs to completion and does not produce spurious calls in normal samples.

## How to Run

<pre>
  docker build -t topmed-mca-loy . 
  docker run -it -v "${PWD}:/workspace" topmed-mca-loy /bin/bash 
  cd .. 
  bash data/scripts/testing/set_up_and_test.sh
</pre>
