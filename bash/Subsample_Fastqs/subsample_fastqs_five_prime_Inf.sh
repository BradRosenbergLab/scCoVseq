#!/bin/bash

#### Count 10X data for scRNA-Seq analysis of SARS CoV-2 infection
#### Sequenced using conventional and modified 5' seq parameters
## Date: 2021-03-01
## Phil Cohen

# Set up path info
PROJECTDIR=#PATH_TO_PROJECT_DIRECTORY
LOGDIR=${PROJECTDIR}/logs

cd ${PROJECTDIR}

# Set target reads per cell to 50,000 reads/cell
TARGET_READS_PER_CELL=50000

################################################################################
# Subsample Conventional 5' Infected Fastqs
################################################################################

# Paths for FASTQ files
CONV_FASTQ_DIR=${PROJECTDIR}/#PATH_TO_CONV5P_FASTQ_DIRECTORY
CONV_INF_FASTQS=${CONV_FASTQ_DIR}/SARSCoV2_5P_Inf_GEX/
WORKING_FASTQ_DIR=${PROJECTDIR}/Subsampled_Fastqs/SARSCoV2_5P_Inf_GEX

# Make working fastq directory if needed
if [[ ! -d $WORKING_FASTQ_DIR ]]; then
  mkdir $WORKING_FASTQ_DIR
fi

###############################################################################
# Combine Fastqs
###############################################################################
COMBINED_FASTQ_DIR=${WORKING_FASTQ_DIR}/Combined

if [[ ! -d $COMBINED_FASTQ_DIR ]]; then
  mkdir $COMBINED_FASTQ_DIR
fi

# R1
# Point to input fastqs
R1_FASTQS="${CONV_INF_FASTQS}/SARSCoV2_5P_Inf_GEX_S1_L001_R1_001.fastq.gz ${CONV_INF_FASTQS}/SARSCoV2_5P_Inf_GEX_S1_L002_R1_001.fastq.gz ${CONV_INF_FASTQS}/SARSCoV2_5P_Inf_GEX_S1_L003_R1_001.fastq.gz ${CONV_INF_FASTQS}/SARSCoV2_5P_Inf_GEX_S1_L004_R1_001.fastq.gz"
# Make name for combined fastqs
COMBINED_R1_FASTQ=${COMBINED_FASTQ_DIR}/SARSCoV2_5P_Inf_GEX_S1_R1_001.fastq.gz
cat ${R1_FASTQS} > ${COMBINED_R1_FASTQ}

# R2
# Point to input fastqs
R2_FASTQS="${CONV_INF_FASTQS}/SARSCoV2_5P_Inf_GEX_S1_L001_R2_001.fastq.gz ${CONV_INF_FASTQS}/SARSCoV2_5P_Inf_GEX_S1_L002_R2_001.fastq.gz ${CONV_INF_FASTQS}/SARSCoV2_5P_Inf_GEX_S1_L003_R2_001.fastq.gz ${CONV_INF_FASTQS}/SARSCoV2_5P_Inf_GEX_S1_L004_R2_001.fastq.gz"
# Make name for combined fastqs
COMBINED_R2_FASTQ=${COMBINED_FASTQ_DIR}/SARSCoV2_5P_Inf_GEX_S1_R2_001.fastq.gz
# Combine fastqs
cat ${R2_FASTQS} > ${COMBINED_R2_FASTQ}

# R2
# Point to input fastqs
I1_FASTQS="${CONV_INF_FASTQS}/SARSCoV2_5P_Inf_GEX_S1_L001_I1_001.fastq.gz ${CONV_INF_FASTQS}/SARSCoV2_5P_Inf_GEX_S1_L002_I1_001.fastq.gz ${CONV_INF_FASTQS}/SARSCoV2_5P_Inf_GEX_S1_L003_I1_001.fastq.gz ${CONV_INF_FASTQS}/SARSCoV2_5P_Inf_GEX_S1_L004_I1_001.fastq.gz"
# Make name for combined fastqs
COMBINED_I1_FASTQ=${COMBINED_FASTQ_DIR}/SARSCoV2_5P_Inf_GEX_S1_I1_001.fastq.gz
# Combine fastqs
cat ${I1_FASTQS} > ${COMBINED_I1_FASTQ}

###############################################################################
# Subsample FASTQs
###############################################################################
SUBSAMPLED_FASTQ_DIR=${WORKING_FASTQ_DIR}/Subsampled

if [[ ! -d $SUBSAMPLED_FASTQ_DIR ]]; then
  mkdir $SUBSAMPLED_FASTQ_DIR
fi

# Calculate Subsampling Depth
# Paths to whitelist
WHITELIST=${PROJECTDIR}/analysis/data/derived_data/Conventional_Inf_Cells.txt
# Determine number of cells
nCells=`wc -l $WHITELIST | awk '{print $1}'`
# Calculate needed read depth
FINAL_DEPTH=`echo "${TARGET_READS_PER_CELL}*${nCells}" | bc -l`

# R1
# Make name for subsampled fastqs
R1_FASTQ_SUBSAMPLED=${SUBSAMPLED_FASTQ_DIR}/Subsampled_SARSCoV2_5P_Inf_GEX_S1_R1_001.fastq.gz
# Subsample combined fastq
seqtk sample -s 123 ${COMBINED_R1_FASTQ} ${FINAL_DEPTH} | gzip > ${R1_FASTQ_SUBSAMPLED}

# R2
# Make name for subsampled fastqs
R2_FASTQ_SUBSAMPLED=${SUBSAMPLED_FASTQ_DIR}/Subsampled_SARSCoV2_5P_Inf_GEX_S1_R2_001.fastq.gz
# Subsample combined fastq
seqtk sample -s 123 ${COMBINED_R2_FASTQ} ${FINAL_DEPTH} | gzip > ${R2_FASTQ_SUBSAMPLED}

# I1
# Make name for subsampled fastqs
I1_FASTQ_SUBSAMPLED=${SUBSAMPLED_FASTQ_DIR}/Subsampled_SARSCoV2_5P_Inf_GEX_S1_I1_001.fastq.gz
# Subsample combined fastq
seqtk sample -s 123 ${COMBINED_I1_FASTQ} ${FINAL_DEPTH} | gzip > ${I1_FASTQ_SUBSAMPLED}
