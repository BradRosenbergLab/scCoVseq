#!/bin/bash

###Phil Cohen
###June 22, 2020
###This script takes as input a viral genome reference and concatenates it with
###a host genome reference. This concatenated genome is then supplied to
###cellranger mkref to generate a reference compatible with downstream
###cellranger functions.

### This function requires cellranger 3.1.0

###Get Viral Reference files and paths
VIRAL_FASTA=$1
VIRAL_GTF=$2
VIRAL_REF_DIR=`cut -d/ -f9 <<< "${VIRAL_FASTA}"`
VIRAL_REF_DIR=`echo $VIRAL_REF_DIR | awk 'BEGIN{FS="."}{print $1}'`

###Get Host Reference files and paths
HOST_REF_DIR=#PATH_TO_YOUR_HOST_REFERENCE
HOST_FASTA=$HOST_REF_DIR/ChlSab1/Chlorocebus_sabaeus.ChlSab1.1.dna.toplevel.fa
HOST_GTF=$HOST_REF_DIR/ChlSab1/Chlorocebus_sabaeus.ChlSab1.1.100.filtered.gtf

###Make output directory
OUTPUT_DIR=#PATH_TO_REFERENCES/$VIRAL_REF_DIR
mkdir $OUTPUT_DIR

###Make combined fasta and gtf
COMBINED_FASTA=`echo "${OUTPUT_DIR}/${VIRAL_REF_DIR}_ChlSab1.fa"`
touch $COMBINED_FASTA
COMBINED_GTF=`echo "${OUTPUT_DIR}/${VIRAL_REF_DIR}_ChlSab1.gtf"`
touch $COMBINED_GTF
FILTERED_GTF=`echo "${OUTPUT_DIR}/${VIRAL_REF_DIR}_ChlSab1_Filtered.gtf"`
touch $FILTERED_GTF

###Combined fasta files
cat $HOST_FASTA > $COMBINED_FASTA
cat $VIRAL_FASTA >> $COMBINED_FASTA

###Combine GTF files
cat $HOST_GTF > $COMBINED_GTF
cat $VIRAL_GTF >> $COMBINED_GTF

###Move to Output Directory
cd $OUTPUT_DIR

###Filter GTF
cellranger mkgtf $COMBINED_GTF $FILTERED_GTF \
--attribute=gene_biotype:protein_coding \
--attribute=gene_biotype:lincRNA \
--attribute=gene_biotype:antisense \

###Make references
cellranger mkref --genome=$VIRAL_REF_DIR \
--fasta=$COMBINED_FASTA \
--genes=$FILTERED_GTF \
--nthreads=8
