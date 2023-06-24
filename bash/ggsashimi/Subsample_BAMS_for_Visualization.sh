#!/bin/bash

#Objective:
# This script will subsample CellRanger output BAMs to 1 x 10^6 viral reads

# Requires samtools v1.11

###Detect input
INPUT_BAM_PATH=$1

###Input depth to subsample reads to 1 x 10^6 Reads
SUBSAMPLED_DEPTH=1000000

###Prepare output paths
SUBSAMPLED_VIRAL_BAM=`echo $INPUT_BAM_PATH | awk 'BEGIN{FS="/"}{print $7 "Downsampled_Viral.bam"}'`
SUBSAMPLED_VIRAL_BAM_DIR=#INPUT_VIRAL_BAM_DIRECTORY
SUBSAMPLED_VIRAL_BAM_PATH=$SUBSAMPLED_VIRAL_BAM_DIR/$SUBSAMPLED_VIRAL_BAM

###Create Subsampled BAM directory if it doesn't already exist
if [[ ! -d $SUBSAMPLED_VIRAL_BAM_DIR ]]; then
  mkdir $SUBSAMPLED_VIRAL_BAM_DIR
fi

###Calculate the number of viral reads in the input file
nReads=`samtools view -c $INPUT_BAM_PATH MN985325.1`

###Calculated fraction needed to get reads to subsampled depth
frac=`echo "${SUBSAMPLED_DEPTH}/${nReads}" | bc -l`

###Subsample input BAM to appropripate depth and showing only viral reads
samtools view -b -s `echo "${frac}"` $INPUT_BAM_PATH MN985325.1 > $SUBSAMPLED_VIRAL_BAM_PATH

###Index the BAM
samtools index $SUBSAMPLED_VIRAL_BAM_PATH
