#!/bin/bash

#### Count 10X data for scRNA-Seq analysis of SARS CoV-2 infection
#### Sequenced using conventional and modified 5' seq parameters
## Date: 2020-07-13
## Phil Cohen

### This script requires cellranger 3.1.0

# Set up path info
PROJECTDIR=#PATH_TO_YOUR_PROJECT
REFDIR=${PROJECTDIR}/analysis/data/raw_data/Combined_References/combined_viral_host_cellranger_indices/Kim_Single_Chrom
COUNTSCRIPTDIR=${PROJECTDIR}/bash/CellRanger_Pipeline/count
LOGDIR=${PROJECTDIR}/CellRanger_Count_Logs

cd ${PROJECTDIR}

# Paths for FASTQ files
FASTQ_DIR=${PROJECTDIR}/Subsampled_Fastqs

THREEP_MOCK_FASTQS=${FASTQ_DIR}/SARSCoV2_3P_Mock_GEX/Subsampled
THREEP_INF_FASTQS=${FASTQ_DIR}/SARSCoV2_3P_Inf_GEX/Subsampled

CONV_MOCK_FASTQS=${FASTQ_DIR}/SARSCoV2_5P_Mock_GEX/Subsampled
CONV_INF_FASTQS=${FASTQ_DIR}/SARSCoV2_5P_Inf_GEX/Subsampled

MOD_MOCK_FASTQS=${FASTQ_DIR}/SARSCoV2_5P_Extended_R1_Mock_GEX/Subsampled
MOD_INF_FASTQS=${FASTQ_DIR}/SARSCoV2_5P_Extended_R1_Infected_GEX/Subsampled

LSF_TEMPLATE=${PROJECTDIR}/bash/CellRanger_Pipeline/lsf.template

###ThreeP Conventional Mock
cellranger count --id `echo ThreeP_Mock_${REFERENCE}` \
		--fastqs=${THREEP_MOCK_FASTQS} \
		--sample=Subsampled_SARSCoV2_3P_Mock_GEX \
		--jobmode=${LSF_TEMPLATE} \
		--transcriptome=${REFDIR} \
		--expect-cells=4000 \
		--nosecondary \
		--maxjobs=200 \
		--chemistry='threeprime'

###ThreeP Conventional Inf
cellranger count --id `echo ThreeP_Inf_${REFERENCE}` \
		--fastqs=${THREEP_INF_FASTQS} \
		--sample=Subsampled_SARSCoV2_3P_Inf_GEX \
		--jobmode=${LSF_TEMPLATE} \
		--transcriptome=${REFDIR} \
		--expect-cells=4000 \
		--nosecondary \
		--maxjobs=200 \
		--chemistry='threeprime'

###FiveP Conventional Mock
cellranger count --id `echo Conventional_Mock_${REFERENCE}` \
		--fastqs=${CONV_MOCK_FASTQS} \
		--sample=Subsampled_SARSCoV2_5P_Mock_GEX \
		--jobmode=${LSF_TEMPLATE} \
		--transcriptome=${REFDIR} \
		--expect-cells=4000 \
		--nosecondary \
		--maxjobs=200 \
		--chemistry='fiveprime'

###FiveP Conventional Inf
cellranger count --id `echo Conventional_Inf_${REFERENCE}` \
		--fastqs=${CONV_INF_FASTQS} \
		--sample=Subsampled_SARSCoV2_5P_Inf_GEX \
		--jobmode=${LSF_TEMPLATE} \
		--transcriptome=${REFDIR} \
		--expect-cells=4000 \
		--nosecondary \
		--maxjobs=200 \
		--chemistry='fiveprime'

##FiveP Extended R1 Mock
cellranger count --id `echo Modified_Mock_${REFERENCE}` \
		--fastqs=${MOD_MOCK_FASTQS} \
		--sample=Subsampled_SARSCoV2_5P_Mock_GEX_reformat \
		--jobmode=${LSF_TEMPLATE} \
		--transcriptome=${REFDIR} \
		--expect-cells=4000 \
		--nosecondary \
		--maxjobs=200 \
		--chemistry='fiveprime'

##FiveP Extended R1 Infected
cellranger count --id `echo Modified_Inf_${REFERENCE}` \
		--fastqs=${MOD_INF_FASTQS} \
		--sample=Subsampled_SARSCoV2_5P_Infected_GEX_reformat \
		--jobmode=${LSF_TEMPLATE} \
		--transcriptome=${REFDIR} \
		--expect-cells=4000 \
		--nosecondary \
		--maxjobs=200 \
		--chemistry='fiveprime'
