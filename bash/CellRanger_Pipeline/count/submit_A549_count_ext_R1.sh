#!/bin/bash

#### Count 10X data for scRNA-Seq analysis of SARS CoV-2 infection (WT, delta ORF6, and ORF6 M58R mutant viruses) of Ace2-expressing A549 cells
#### Sequenced using Extended R1 method
## Date: 2022-07-04
## Phil Cohen

### This script requires cellranger 3.1.0

# Set up path info
export PROJECTDIR=#PATH_TO_YOUR_PROJECT
export REFDIR=${PROJECTDIR}/analysis/data/raw_data/Combined_References/combined_viral_host_cellranger_indices/Human_CoV2_Reference/hg38_Kim_Single_Chrom
export COUNTSCRIPTDIR=${PROJECTDIR}/bash/CellRanger_Pipeline/count
export LOGDIR=${PROJECTDIR}/CellRanger_Count_Logs
export LSF_TEMPLATE=${PROJECTDIR}/bash/CellRanger_Pipeline/lsf.template

cd ${PROJECTDIR}

# Paths for FASTQ files for Run 1
export FASTQ_DIR=#PATH_TO_YOUR_FASTQS
export A549_MOCK_1_FASTQS=${FASTQ_DIR}/M1
export A549_M58R_FASTQS=${FASTQ_DIR}/M58R

# Submit count scripts
## Mock_1
cellranger count --id A549_Mock_1 \
		--fastqs=${A549_MOCK_1_FASTQS} \
		--sample=M1_reformat \
		--jobmode=${LSF_TEMPLATE} \
		--transcriptome=${REFDIR} \
		--expect-cells=4000 \
		--nosecondary \
		--maxjobs=200 \
		--chemistry='fiveprime'

## M58R infection
cellranger count --id count_A549_M58R \
		--fastqs=${A549_M58R_FASTQS} \
		--sample=M58R_reformat \
		--jobmode=${LSF_TEMPLATE} \
		--transcriptome=${REFDIR} \
		--expect-cells=4000 \
		--nosecondary \
		--maxjobs=200 \
		--chemistry='fiveprime'
