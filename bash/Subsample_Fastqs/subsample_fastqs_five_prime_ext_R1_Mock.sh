#!/bin/bash

#### Count 10X data for scRNA-Seq analysis of SARS CoV-2 infection
#### Sequenced using conventional and modified 5' seq parameters
## Date: 2021-03-01
## Phil Cohen

module load seqtk/1.2

# Set up path info
PROJECTDIR=#PATH_TO_PROJECT_DIRECTORY
LOGDIR=${PROJECTDIR}/logs

cd ${PROJECTDIR}

# Set target reads per cell to 50,000 reads/cell
TARGET_READS_PER_CELL=50000

################################################################################
# Subsample Modified 5' Mock Fastqs
################################################################################

# Paths for FASTQ files
MOD_FASTQ_DIR=${PROJECTDIR}/SARSCoV2_Sequencing/Reformatted/
MOD_MOCK_FASTQS=${MOD_FASTQ_DIR}/SARSCoV2_5P_Extended_R1_Mock_GEX
WORKING_FASTQ_DIR=${PROJECTDIR}/Subsampled_Fastqs/SARSCoV2_5P_Extended_R1_Mock_GEX

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
R1_FASTQS="${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S1_L001_R1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S1_L002_R1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S1_L003_R1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S1_L004_R1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S2_L001_R1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S2_L002_R1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S2_L003_R1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S2_L004_R1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S3_L001_R1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S3_L002_R1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S3_L003_R1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S3_L004_R1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S4_L001_R1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S4_L002_R1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S4_L003_R1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S4_L004_R1_001.fastq.gz"
# Make name for combined fastqs
COMBINED_R1_FASTQ=${COMBINED_FASTQ_DIR}/SARSCoV2_5P_Mock_GEX_reformat_R1_001.fastq.gz
cat ${R1_FASTQS} > ${COMBINED_R1_FASTQ}

# R2
# Point to input fastqs
R2_FASTQS="${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S1_L001_R2_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S1_L002_R2_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S1_L003_R2_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S1_L004_R2_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S2_L001_R2_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S2_L002_R2_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S2_L003_R2_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S2_L004_R2_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S3_L001_R2_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S3_L002_R2_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S3_L003_R2_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S3_L004_R2_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S4_L001_R2_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S4_L002_R2_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S4_L003_R2_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S4_L004_R2_001.fastq.gz"
# Make name for combined fastqs
COMBINED_R2_FASTQ=${COMBINED_FASTQ_DIR}/SARSCoV2_5P_Mock_GEX_reformat_R2_001.fastq.gz
# Combine fastqs
cat ${R2_FASTQS} > ${COMBINED_R2_FASTQ}

# R2
# Point to input fastqs
I1_FASTQS="${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S1_L001_I1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S1_L002_I1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S1_L003_I1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S1_L004_I1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S2_L001_I1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S2_L002_I1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S2_L003_I1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S2_L004_I1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S3_L001_I1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S3_L002_I1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S3_L003_I1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S3_L004_I1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S4_L001_I1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S4_L002_I1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S4_L003_I1_001.fastq.gz ${MOD_MOCK_FASTQS}/SARSCoV2_5P_Mock_GEX_reformat_S4_L004_I1_001.fastq.gz"
# Make name for combined fastqs
COMBINED_I1_FASTQ=${COMBINED_FASTQ_DIR}/SARSCoV2_5P_Mock_GEX_reformat_I1_001.fastq.gz
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
WHITELIST=${PROJECTDIR}/analysis/data/derived_data/Modified_Mock_Cells.txt
# Determine number of cells
nCells=`wc -l $WHITELIST | awk '{print $1}'`
# Calculate needed read depth
FINAL_DEPTH=`echo "${TARGET_READS_PER_CELL}*${nCells}" | bc -l`

# R1
# Make name for subsampled fastqs
R1_FASTQ_SUBSAMPLED=${SUBSAMPLED_FASTQ_DIR}/Subsampled_SARSCoV2_5P_Mock_GEX_reformat_R1_001.fastq.gz
# Subsample combined fastq
seqtk sample -s 123 ${COMBINED_R1_FASTQ} ${FINAL_DEPTH} | gzip > ${R1_FASTQ_SUBSAMPLED}

# R2
# Make name for subsampled fastqs
R2_FASTQ_SUBSAMPLED=${SUBSAMPLED_FASTQ_DIR}/Subsampled_SARSCoV2_5P_Mock_GEX_reformat_R2_001.fastq.gz
# Subsample combined fastq
seqtk sample -s 123 ${COMBINED_R2_FASTQ} ${FINAL_DEPTH} | gzip > ${R2_FASTQ_SUBSAMPLED}

# I1
# Make name for subsampled fastqs
I1_FASTQ_SUBSAMPLED=${SUBSAMPLED_FASTQ_DIR}/Subsampled_SARSCoV2_5P_Mock_GEX_reformat_I1_001.fastq.gz
# Subsample combined fastq
seqtk sample -s 123 ${COMBINED_I1_FASTQ} ${FINAL_DEPTH} | gzip > ${I1_FASTQ_SUBSAMPLED}
