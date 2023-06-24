#!/bin/bash

INPUT_FASTA=/sc/arion/projects/brrbrtseq/SARS-CoV-2_5P_10X_Sequencing/analysis/data/raw_data/Combined_References/Modified_CoV2_Reference_Files/Kim_Single_Chrom.fa
INPUT_GTF=/sc/arion/projects/brrbrtseq/SARS-CoV-2_5P_10X_Sequencing/analysis/data/raw_data/Combined_References/Modified_CoV2_Reference_Files/Kim_Single_Chrom.gtf

STAR --runThreadN 6 \
  --runMode genomeGenerate \
  --genomeDir ./SARSCoV2_Genome \
  --genomeSAindexNbases 6 \
  --genomeFastaFiles $INPUT_FASTA \
  --sjdbGTFfile $INPUT_GTF \
  --sjdbOverhang 118 \
  --genomeSAsparseD 3
