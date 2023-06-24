#!/bin/bash

#Objective:
# This script will perform a modified mapping and quantification strategy
# for 10X scRNAseq datasets intended to study SARS-CoV-2, or other coronavirus,
# infection.

# The pipeline can be broken down into three parts:
# 1. Subset reads mapping to the SARS-CoV-2 viral transcriptome reference and
#    filter viral reads that do not contain leader junctions or
#    part of ORF1a or ORF1ab (Region 1:21,555, the final base of the STOP codon of ORF1ab)
# 2. Make a Host BAM that excludes all Viral Reads
# 3. Merge filtered viral BAMs with host BAM
# 4. Quantify UMIs to generate a Gene x Cell Matrix with umi_tools

# Code requires:
## samtools v1.11
## python v3.7.3
## R v3.5.3

##PART 0: Prepare Input Files
###Prepare input files
INPUT_BAM_PATH=$1
OUTPUT_PATH=#YOUR_OUTPUT_PATH
OUTPUT_BAM_PATH=$OUTPUT_PATH/BAMs_from_Quantification_Pipelines

if [[ ! -d $OUTPUT_BAM_PATH ]]; then
  mkdir $OUTPUT_BAM_PATH
fi

###Prepare output files
CELLRANGER_OUTPUT_BAM_DIR=${OUTPUT_BAM_PATH}/CellRanger_Output_BAMs
CELLRANGER_OUTPUT_BAM_FILE=`echo $INPUT_BAM_PATH | awk 'BEGIN{FS="/"}{print $7 ".bam"}'`
CELLRANGER_OUTPUT_BAM_PATH=$CELLRANGER_OUTPUT_BAM_DIR/$CELLRANGER_OUTPUT_BAM_FILE

##Create Cellranger Output BAM directory if it doesn't already exist
if [[ ! -d $CELLRANGER_OUTPUT_BAM_DIR ]]; then
  mkdir $CELLRANGER_OUTPUT_BAM_DIR
fi

###Copy CellRanger Output BAM to CellRanger_Output_BAMs subdirectory
cp $INPUT_BAM_PATH $CELLRANGER_OUTPUT_BAM_PATH

###Index CellRanger Output BAM
samtools index $CELLRANGER_OUTPUT_BAM_PATH

## PART 1: Split Genomic and Subgenomic Viral Reads
### All genomic reads must be contiguous and ungapped and
### map anywhere from the start of the genome to the 3'
### end of ORF1ab

### All subgenomic reads must contain a gap and map in
### part to the leader region. We define the leader region
### as bp 1-80 of the SARS-CoV-2 genome. The longest annotated
### leader region of the Kim sgmRNAs is 76 bp (citation below),
### but the length of the leader region varies between sgmRNAs,
### so we will use a wider filter.
### Kim, D., Lee, J.-Y., Yang, J.-S., Kim, J.W., Kim, V.N., and Chang, H. (2020).
### The architecture of SARS-CoV-2 transcriptome. Cell 181, 914-921.e10.
### 10.1016/j.cell.2020.04.011.


###Create Junctional directory if it doesn't already exist
JUNCTIONAL_DIR=${OUTPUT_BAM_PATH}/Junctional
if [[ ! -d $JUNCTIONAL_DIR ]]; then
  mkdir $JUNCTIONAL_DIR
fi

###Create Viral BAM directory if it doesn't already exist
VIRAL_BAM_DIR=${JUNCTIONAL_DIR}/Viral_BAMS
if [[ ! -d $VIRAL_BAM_DIR ]]; then
  mkdir $VIRAL_BAM_DIR
fi

###Prepare output files
GENOME_VIRAL_BAM=`echo "${CELLRANGER_OUTPUT_BAM_FILE}" | awk 'BEGIN{FS="."} {print $1 "_Viral_Genome.bam"}'`
GENOME_VIRAL_BAM_PATH=$VIRAL_BAM_DIR/$GENOME_VIRAL_BAM

SUBGENOME_VIRAL_BAM=`echo "${CELLRANGER_OUTPUT_BAM_FILE}" | awk 'BEGIN{FS="."} {print $1 "_Viral_Subgenome.bam"}'`
SUBGENOME_VIRAL_BAM_TEMP=`echo "${CELLRANGER_OUTPUT_BAM_FILE}" | awk 'BEGIN{FS="."} {print $1 "_Viral_Subgenome_temp.bam"}'`

SUBGENOME_VIRAL_BAM_PATH=$VIRAL_BAM_DIR/$SUBGENOME_VIRAL_BAM
SUBGENOME_VIRAL_BAM_TEMP_PATH=$VIRAL_BAM_DIR/$SUBGENOME_VIRAL_BAM_TEMP

VIRAL_BAM=`echo "${CELLRANGER_OUTPUT_BAM_FILE}" | awk 'BEGIN{FS="."} {print $1 "_Viral.bam"}'`
VIRAL_BAM_PATH=$VIRAL_BAM_DIR/$VIRAL_BAM

SORTED_VIRAL_BAM=`echo "${CELLRANGER_OUTPUT_BAM_FILE}" | awk 'BEGIN{FS="."} {print $1 "_Viral_sorted.bam"}'`
SORTED_VIRAL_BAM_PATH=$VIRAL_BAM_DIR/$SORTED_VIRAL_BAM

### Filter to Genomic RNA
### Subset to only ungapped viral reads mapping to ORF1a/ORF1ab Region
samtools view -h $CELLRANGER_OUTPUT_BAM_PATH MN985325.1:1-21,549 | awk '($0 ~ /^@/) || (!($6 ~ /N/))' | samtools view -b > $GENOME_VIRAL_BAM_PATH
samtools index $GENOME_VIRAL_BAM_PATH

### Filter to Subgenomic RNA
### Subset to only gapped viral reads that contain leader sequence and sgmRNA sequence
samtools view -h $CELLRANGER_OUTPUT_BAM_PATH MN985325.1:1-80 | awk '($0 ~ /^@/) || ($6 ~ /N/)' | samtools view -b > $SUBGENOME_VIRAL_BAM_TEMP_PATH
samtools index $SUBGENOME_VIRAL_BAM_TEMP_PATH
samtools view -b $SUBGENOME_VIRAL_BAM_TEMP_PATH MN985325.1:21,550-29,882 > $SUBGENOME_VIRAL_BAM_PATH
rm $SUBGENOME_VIRAL_BAM_TEMP_PATH

### Merge genomic and subgenomic viral bams into viral bam
samtools merge -f $VIRAL_BAM_PATH $SUBGENOME_VIRAL_BAM_PATH $GENOME_VIRAL_BAM_PATH
samtools index $VIRAL_BAM_PATH

### Sort bam
samtools sort -o $SORTED_VIRAL_BAM_PATH $VIRAL_BAM_PATH
samtools index $SORTED_VIRAL_BAM_PATH

##PART 2: Remove Reads Aligning to Virus from BAM to make Host BAM
###Prepare output files
HOST_BAM=`echo "${CELLRANGER_OUTPUT_BAM_FILE}" | awk 'BEGIN{FS="."} {print $1 "_Host.bam"}'`
HOST_BAM_DIR=${JUNCTIONAL_DIR}/Host_BAMS
HOST_BAM_PATH=$HOST_BAM_DIR/$HOST_BAM

SORTED_HOST_BAM=`echo "${CELLRANGER_OUTPUT_BAM_FILE}" | awk 'BEGIN{FS="."} {print $1 "_Host_sorted.bam"}'`
SORTED_HOST_BAM_PATH=$HOST_BAM_DIR/$SORTED_HOST_BAM

###Create Host BAM directory if it doesn't already exist
if [[ ! -d $HOST_BAM_DIR ]]; then
  mkdir $HOST_BAM_DIR
fi

###Select Reads that do not align to virus and output to Host BAM
samtools view -b $CELLRANGER_OUTPUT_BAM_PATH -L filter.bed -U $HOST_BAM_PATH -o /dev/null

###Sort BAM
samtools sort -o $SORTED_HOST_BAM_PATH $HOST_BAM_PATH

###Index host BAM
samtools index $SORTED_HOST_BAM_PATH

##PART 3: Merge Viral and Host BAMS
###Prepare output files
MERGED_BAM_DIR=${JUNCTIONAL_DIR}/Merged_Virus_Host_BAMS
MERGED_BAM=`echo "${CELLRANGER_OUTPUT_BAM_FILE}" | awk 'BEGIN{FS="."} {print $1 "_Merged.bam"}'`
MERGED_BAM_PATH=$MERGED_BAM_DIR/$MERGED_BAM

SORTED_MERGED_BAM=`echo "${CELLRANGER_OUTPUT_BAM_FILE}" | awk 'BEGIN{FS="."} {print $1 "_Merged_sorted.bam"}'`
SORTED_MERGED_BAM_PATH=$MERGED_BAM_DIR/$SORTED_MERGED_BAM

###Create Junctional BAM directory if it doesn't already exist
if [[ ! -d $MERGED_BAM_DIR ]]; then
  mkdir $MERGED_BAM_DIR
fi

###Merge Host and Virus BAM
samtools merge -fh $SORTED_HOST_BAM_PATH $MERGED_BAM_PATH $SORTED_HOST_BAM_PATH $SORTED_VIRAL_BAM_PATH

###Sort merged bam
samtools sort -o $SORTED_MERGED_BAM_PATH $MERGED_BAM_PATH

###Index merged bam
samtools index $SORTED_MERGED_BAM_PATH

##PART 4: Quantify UMIs of Junctional BAM
###Prepare output files
GENE_CELL_MATRIX=`echo "${CELLRANGER_OUTPUT_BAM_FILE}" | awk 'BEGIN{FS="."} {print $1 "_Gene_Cell_Matrix.tsv.gz"}'`
GENE_CELL_MATRIX_DIR=${OUTPUT_PATH}/scCoVseq_Gene_x_Cell_Matrices
GENE_CELL_MATRIX_PATH=$GENE_CELL_MATRIX_DIR/$GENE_CELL_MATRIX

###Create Gene Cell Matrix directory if it doesn't already exist
if [[ ! -d $GENE_CELL_MATRIX_DIR ]]; then
  mkdir $GENE_CELL_MATRIX_DIR
fi

###Calculate Gene Counts and output Gene x Cell Matrix
umi_tools count --extract-umi-method=tag --umi-tag=UB --cell-tag=CB --method=unique --per-gene --gene-tag=GN --per-cell --unmapped-reads=discard --wide-format-cell-counts -I $SORTED_MERGED_BAM_PATH -S $GENE_CELL_MATRIX_PATH

###Convert tsv.gz to sparse matrix and save as rds
Rscript ./convert_tsv_to_mtx.R $GENE_CELL_MATRIX_PATH
rm $GENE_CELL_MATRIX_PATH

echo "Success"
