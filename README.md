# wellDR-seq

This repository contains scripts and meta data metrics files used in the wellDR-seq paper (Coalescing single cell genomes and transcriptomes to decode breast cancer progression). wellDR-seq is a high-throughput single cell sequencing technology that enables simultaneous profiling the whole geneome and transcrptome from thousands of single cells.

## Project Overview

This project provides a complete analysis pipeline including data preprocessing, quality control, copy number analysis, gene expression analysis, and multi-omics integration analysis for wellDR-seq data.

## Directory Structure

### scripts/
R scripts directory contains codes that can reproduce the analysis and figures from the wellDR-seq paper:

- **0_wdr_functions.R**: Core function library containing main functions for wellDR-seq analysis
- **DNA_QC.R**: DNA data quality control and analysis script
- **RNA_QC.R**: RNA data quality control and analysis script
- **gene_dosage.R**: Gene dosage analysis script
- **sporadic_CNA.R**: Sporadic copy number aberrations analysis
- **resolution_corr_downsampling.R**: Resolution correlation and downsampling analysis
- **mutation_calling.sh**: Mutation calling shell script
- **P1-P12_*_process.R**: Processing scripts for 12 different samples
- **wdr_mda231_process.R**: MDA231 cell line wellDR-seq data processing
- **rna_fq_preprocess/**: RNA sequencing data preprocessing scripts
  - **scDR-fq-convert-v4.pl**: Sequencing data format conversion script
  - **star-smart.sh**: STAR alignment and SMART-seq processing script
  - **com_index_wafer.txt**: Chip index file

### pre_load_data/
Data files required for R script execution:

- **pre_load_data.rda**: Pre-loaded R data objects
- **bin_coords.csv**: Genomic bin coordinate information
- **hg19_gene_binpos_map.tsv**: hg19 genome gene position mapping
- **chr_*.txt**: Chromosome-related information files
- **allcell_markers.csv**: Cell marker gene information
- **HBCA_Epi_markers.csv**: HBCA epithelial cell marker genes
- **wafer_match_list.csv**: Chip matching list
- **bin.boundaries.50k.bowtie.k50.sorted.bin_not_removed.txt**: 50kb bin boundary file
- **50kb_excluded.txt**: Excluded regions for 50kb bin resolution


## Data Source

The original sequencing data from this study has been deposited to the NCBI SRA: PRJNA1086561 and the NCBI GEO: GSE261713.

## Citation

If you use the wellDR-seq method in your research, please cite the relevant paper.

## Contact

For any additional information, please email the corresponding authors.

