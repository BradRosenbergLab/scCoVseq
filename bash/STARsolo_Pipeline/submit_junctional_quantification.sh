#!/bin/bash

###Store genomes to test
INPUT_REFERENCE=/sc/arion/projects/brrbrtseq/SARS-CoV-2_5P_10X_Sequencing/bash/Junctional_Quantification/SARSCoV2_Genome
INPUT_BAM_PATH=/sc/arion/projects/brrbrtseq/SARS-CoV-2_5P_10X_Sequencing/Modified_Inf_Kim_Single_Chrom/outs/possorted_genome_bam.bam
INPUT_WHITELIST=/sc/arion/projects/brrbrtseq/SARS-CoV-2_5P_10X_Sequencing/analysis/data/derived_data/Modified_Inf_Cells.txt

STAR --runThreadN 6 \
    --genomeDir $INPUT_REFERENCE \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist None \
    --soloFeatures SJ \
    --readFilesIn $INPUT_BAM_PATH \
    --readFilesType SAM SE \
    --readFilesCommand samtools view -F 0x100 \
    --soloInputSAMattrBarcodeSeq CR UR \
    --soloInputSAMattrBarcodeQual CY UY

echo "Alignment complete"
