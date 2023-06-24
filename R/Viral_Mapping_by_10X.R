library(paletteer)
library(mgsub)
library(Rsamtools)
library(tidyverse)
library(future.apply)

### Specify parameters for BAM importing
BAM_Params <- ScanBamParam(
  tag = c("GN", "NH"),
  which = GRanges(
    seqnames = "MN985325.1",
    ranges = IRanges(
      start = 1,
      end = 29882
    )
  ),
  what = c("qname", "pos", "cigar")
)

### Specify paths to bams
project_path <- #PATH_TO_YOUR_PROJECT
viral_bam_path <- paste0(project_path, "/analysis/data/raw_data/BAMs_from_Quantification_Pipelines/Junctional/Viral_BAMS")
cellranger_bam_path <- paste0(project_path, "/analysis/data/raw_data/BAMs_from_Quantification_Pipelines/CellRanger_Output_BAMs")

## 10X 3' BAMs
### Specify specific bams
input_BAM <- paste0(cellranger_bam_path, "/ThreeP_Inf_Kim_Single_Chrom.bam")
subgenome_BAM <- paste0(viral_bam_path, "/ThreeP_Inf_Kim_Single_Chrom_Viral_Subgenome.bam")
genome_BAM <- paste0(viral_bam_path, "/ThreeP_Inf_Kim_Single_Chrom_Viral_Genome.bam")

libraries <- c("10X 3'", "10X 5'", "10X 5' Extended R1")
plan(multisession, workers = 4)
Viral_Gene_Multimapping <- future_lapply(libraries, function(library) {
  if (library == "10X 3'") {
    input_BAM <- paste0(cellranger_bam_path, "/ThreeP_Inf_Kim_Single_Chrom.bam")
    subgenome_BAM <- paste0(viral_bam_path, "/ThreeP_Inf_Kim_Single_Chrom_Viral_Subgenome.bam")
    genome_BAM <- paste0(viral_bam_path, "/ThreeP_Inf_Kim_Single_Chrom_Viral_Genome.bam")
  } else if (library == "10X 5'") {
    input_BAM <- paste0(cellranger_bam_path, "/Conventional_Inf_Kim_Single_Chrom.bam")
    subgenome_BAM <- paste0(viral_bam_path, "/Conventional_Inf_Kim_Single_Chrom_Viral_Subgenome.bam")
    genome_BAM <- paste0(viral_bam_path, "/Conventional_Inf_Kim_Single_Chrom_Viral_Genome.bam")
  } else if (library == "10X 5' Extended R1") {
    input_BAM <- paste0(cellranger_bam_path, "/Modified_Inf_Kim_Single_Chrom.bam")
    subgenome_BAM <- paste0(viral_bam_path, "/Modified_Inf_Kim_Single_Chrom_Viral_Subgenome.bam")
    genome_BAM <- paste0(viral_bam_path, "/Modified_Inf_Kim_Single_Chrom_Viral_Genome.bam")
  }
  ### Import whole genome bam
  whole_genome_bam <- scanBam(
    file = input_BAM,
    param = BAM_Params
  )

  ### Convert bam to dataframe and set read_type to ambiguous
  whole_genome_bam_df <- data.frame(
    Read_Name = whole_genome_bam$`MN985325.1:1-29882`$qname,
    Start_Site = whole_genome_bam$`MN985325.1:1-29882`$pos,
    CIGAR = whole_genome_bam$`MN985325.1:1-29882`$cigar,
    Gene = whole_genome_bam$`MN985325.1:1-29882`$tag$GN,
    Number_of_Maps = whole_genome_bam$`MN985325.1:1-29882`$tag$NH,
    Read_Type = "Ambiguous"
  )


  ### Subgenome
  subgenome_BAM <- scanBam(
    file = subgenome_BAM,
    param = BAM_Params
  )

  ### Convert bam to dataframe and set read_type to Leader-sgmRNA Junction-Spanning
  subgenome_BAM_df <- data.frame(
    Read_Name = subgenome_BAM$`MN985325.1:1-29882`$qname,
    Start_Site = subgenome_BAM$`MN985325.1:1-29882`$pos,
    CIGAR = subgenome_BAM$`MN985325.1:1-29882`$cigar,
    Gene = subgenome_BAM$`MN985325.1:1-29882`$tag$GN,
    Number_of_Maps = subgenome_BAM$`MN985325.1:1-29882`$tag$NH,
    Read_Type = "Leader-sgmRNA\nJunction-Spanning"
  )

  ### Genome
  genome_BAM <- scanBam(
    file = genome_BAM,
    param = BAM_Params
  )

  ### Convert bam to dataframe and set read_type to Contiguous ORF1a/b (gRNA)
  genome_BAM_df <- data.frame(
    Read_Name = genome_BAM$`MN985325.1:1-29882`$qname,
    Start_Site = genome_BAM$`MN985325.1:1-29882`$pos,
    CIGAR = genome_BAM$`MN985325.1:1-29882`$cigar,
    Gene = genome_BAM$`MN985325.1:1-29882`$tag$GN,
    Number_of_Maps = genome_BAM$`MN985325.1:1-29882`$tag$NH,
    Read_Type = "Contiguous ORF1a/b\n(gRNA)"
  )

  ### Combine read names for unambiguous reads
  quantified_reads <- union(genome_BAM_df$Read_Name, subgenome_BAM_df$Read_Name)

  ### Select reads that are not unambiguous
  excluded_reads <- whole_genome_bam_df$Read_Name[!whole_genome_bam_df$Read_Name %in% quantified_reads]

  ### Subset whole genome bam to excluded reads
  excluded_reads_BAM_df <- subset(whole_genome_bam_df, Read_Name %in% excluded_reads)

  ### Combine into one dataframe
  output <- rbind(
    subgenome_BAM_df,
    genome_BAM_df,
    excluded_reads_BAM_df
  )

  ### Set sequencing method
  output$Sequencing_Method <- library

  return(output)
})

### Combine output into one dataframe
Viral_Gene_Multimapping <- bind_rows(Viral_Gene_Multimapping)

### Calculate the number of genes each read mapped to
Viral_Gene_Multimapping$Number_of_Genes <- sapply(str_split(Viral_Gene_Multimapping$Gene, pattern = ";"), length)

### Set the levels for read types
Viral_Gene_Multimapping$Read_Type <- factor(Viral_Gene_Multimapping$Read_Type,
  levels = c(
    "Ambiguous",
    "Contiguous ORF1a/b\n(gRNA)",
    "Leader-sgmRNA\nJunction-Spanning"
  )
)

### Reformat data
Viral_Gene_Multimapping <- Viral_Gene_Multimapping %>%
  ### Calculate number of reads per read type per number of genes per library
  group_by(Read_Type, Sequencing_Method, Number_of_Genes) %>%
  add_tally(name = "nReads") %>%
  ### Calculate total number of reads in a library
  group_by(Sequencing_Method) %>%
  add_tally(name = "Reads_per_Library") %>%
  ### Remove unneeded columns
  select(Read_Type, Sequencing_Method, Number_of_Genes, nReads, Reads_per_Library) %>%
  ### Select only distinct rows
  distinct() %>%
  ### Calculate the fraction of reads mapping to a single gene by read type per library
  mutate(Fraction_of_Reads = nReads / Reads_per_Library) %>%
  ### Calculate the reads per million mapping to a single gene by read type per library
  mutate(Reads_per_million = Fraction_of_Reads * 10^6)

output_path <- paste0(project_path, "/analysis/figures/Vectorized/Figure_1")

saveRDS(Viral_Gene_Multimapping, file = paste0(output_path, "/Viral_Gene_Multimapping.rds"))
