#!/bin/bash

### This script submits a periscope job for each sample in the SAMPLES variable
### Phil Cohen, PhD
### Last Edited: 8/30/2022

### Prepare a list of samples
export SAMPLES="Mock_2_S1 WT_1_S2 WT_2_S3 dORF6_1_S4 dORF6_2_S5 dORF6_3_S6 M58R_1_S7 M58R_2_S8 M58R_3_S9 Mock_1_S10 Mock_3_S11 WT_3_S12"

### Prepare paths
export PROJECTDIR=#PATH_TO_YOUR_PROJECT
export FASTQ_1_PATH=${PROJECTDIR}/HW5KWAFX3
export FASTQ_2_PATH=${PROJECTDIR}/HW5T3AFX3
export OUTPUT_PATH=${PROJECTDIR}/analysis/data/raw_data/Bulk_RNAseq/Periscope
export COMBINED_FASTQ_PATH=${PROJECTDIR}/analysis/data/raw_data/Bulk_RNAseq/Combined_fastqs

### Make an output directory if it doesn't already exist
if [[ ! -d ${OUTPUT_PATH} ]]; then
  mkdir ${OUTPUT_PATH}
fi

### Make a combined fastq directory if it doesn't already exist
if [[ ! -d ${COMBINED_FASTQ_PATH} ]]; then
  mkdir ${COMBINED_FASTQ_PATH}
fi

### Loop through samples and submit an alignment job per sample
for i in $SAMPLES; do
  ### Make an output directory per sample if it doesn't already exist
  if [[ ! -d ${OUTPUT_PATH}/${i} ]]; then
    mkdir ${OUTPUT_PATH}/${i}
  fi

  ### Combine R1 fastqs if it doesn't already exist
  if [[ ! -f ${COMBINED_FASTQ_PATH}/${i}_R1_001.fastq.gz ]]; then
    cat ${FASTQ_1_PATH}/${i}_R1_001.fastq.gz ${FASTQ_2_PATH}/${i}_R1_001.fastq.gz > ${COMBINED_FASTQ_PATH}/${i}_R1_001.fastq.gz
  fi

  ### Combine R2 fastqs if it doesn't already exist
  if [[ ! -f ${COMBINED_FASTQ_PATH}/${i}_R2_001.fastq.gz ]]; then
    cat ${FASTQ_1_PATH}/${i}_R2_001.fastq.gz ${FASTQ_2_PATH}/${i}_R2_001.fastq.gz > ${COMBINED_FASTQ_PATH}/${i}_R2_001.fastq.gz
  fi

  ### Submit an alignment per sample
  cd ${OUTPUT_PATH}/${i};
  periscope \
      --fastq ${COMBINED_FASTQ_PATH}/${i}_R1_001.fastq.gz ${COMBINED_FASTQ_PATH}/${i}_R2_001.fastq.gz \
      --output-prefix ${OUTPUT_PATH}/${i}/${i}_periscope \
      --sample ${i} \
      --tmp ${OUTPUT_PATH}/${i}/ \
      --technology illumina \
      --threads 16
done
echo "Success"
