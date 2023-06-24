#!/bin/bash

#### Count 10X data for scRNA-Seq analysis of SARS CoV-2 infection
#### Sequenced using conventional and modified 5' seq parameters
## Date: 2020-07-13
## Phil Cohen

### This script requires cellranger 3.1.0

# Set up path info
export PROJECTDIR=#PATH_TO_YOUR_PROJECT
export REFERENCE=${PROJECTDIR}/analysis/data/raw_data/Combined_References/combined_viral_host_cellranger_indices/NCBI_Single_Chrom
export COUNTSCRIPTDIR=${PROJECTDIR}/bash/CellRanger_Pipeline/count
export LOGDIR=${PROJECTDIR}/CellRanger_Count_Logs

cd ${PROJECTDIR}

# Paths for FASTQ files
export CONV_FASTQS=${PROJECTDIR}/HVH25BGXF/outs/fastq_path/HVH25BGXF/
export CONV_FASTQS_2=${PROJECTDIR}/HMM5FAFX2/outs/fastq_path/HMM5FAFX2/
export MOD_GEX=${PROJECTDIR}/SARSCoV2_Sequencing/Reformatted/
export LSF_TEMPLATE=${PROJECTDIR}/bash/CellRanger_Pipeline/lsf.template

##Modified Mock
cellranger count --id `echo Modified_Mock_${REFERENCE}` \
		--fastqs=${MOD_GEX} \
		--sample=SARSCoV2_5P_Mock_GEX_reformat\
		--jobmode=${LSF_TEMPLATE} \
		--transcriptome=${REFERENCE} \
		--expect-cells=4000 \
		--nosecondary \
		--maxjobs=200 \
		--chemistry='fiveprime'

##Modified Infected
cellranger count --id `echo Modified_Inf_${REFERENCE}` \
		--fastqs=${MOD_GEX} \
		--sample=SARSCoV2_5P_Infected_GEX_reformat \
		--jobmode=${LSF_TEMPLATE} \
		--transcriptome=${REFERENCE} \
		--expect-cells=4000 \
		--nosecondary \
		--maxjobs=200 \
		--chemistry='fiveprime'

###Conventional Mock
cellranger count --id `echo Conventional_Mock_${REFERENCE}` \
		--fastqs=${CONV_FASTQS},${CONV_FASTQS_2} \
		--sample=SARSCoV2_5P_Mock_GEX \
		--jobmode=${LSF_TEMPLATE} \
		--transcriptome=${REFERENCE} \
		--expect-cells=4000 \
		--nosecondary \
		--maxjobs=200 \
		--chemistry='fiveprime'

###Conventional Inf
cellranger count --id `echo Conventional_Inf_${REFERENCE}` \
		--fastqs=${CONV_FASTQS},${CONV_FASTQS_2} \
		--sample=SARSCoV2_5P_Inf_GEX \
		--jobmode=${LSF_TEMPLATE} \
		--transcriptome=${REFERENCE} \
		--expect-cells=4000 \
		--nosecondary \
		--maxjobs=200 \
		--chemistry='fiveprime'
