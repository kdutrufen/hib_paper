library(ballgown)
library(CEMiTool)
library(clusterProfiler)
library(cowplot)
library(data.table)
library(DESeq2)
library(EDASeq)
library(edgeR)
library(enrichR)
library(eply)
library(formattable)
library(ggplot2)
library(gplots)
library(ggpubr)
library(ggrepel)
library(matrixStats)
library(patchwork)
library(pathview)
library(ReactomePA)
library(RColorBrewer)
library(RUVSeq)
library(STRINGdb)
library(stringr)
library(tidyverse)
library(VennDiagram)
library(vroom)
library(xtable)

source_path <- "C:/Users/Carlos/Dropbox/"
# source_path <- "~/Dropbox/"

source(paste0(source_path, "kdutrufenr/kdutrufenr.R"))

path_to_files <- paste0(source_path, "Hib_Data/")

# Protein list
Hib_proteins <- path_to_files %>%
  paste0("Haemophilus_influenzae_10810.GCF_000210875.1.proteintable.txt") %>%
  read_delim(delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  as.data.frame()

# Read in Hib Matrix form Hisat alignment
# Count matrix v3 (Hisat k -5, featureCounts -M --primary)
Hib_count_data <- path_to_files %>%
  paste0("Hib_Hisat2_featureCounts_primary_table.txt") %>%
  read_delim("\t", escape_double = FALSE, trim_ws = TRUE, skip = 1) %>%
  as.data.frame()
Hib_count_data <- Hib_count_data %>%
  dplyr::filter(!str_detect(Geneid, pattern = "rrf")) %>%
  column_to_rownames("Geneid") %>%
  dplyr::select(-c(Chr, Start, End, Strand, Length)) %>%
  setNames(sort(as.vector(outer(paste0("S0", 1:6), paste0("B0", 1:4), paste0))))

Hib_count_data <- Hib_count_data %>% dplyr::select(-(Hib_count_data %>% colnames() %>% str_which(pattern = "B01")))

samples_to_remove <- c("S03B02")
Hib_count_data <- Hib_count_data %>% dplyr::select(-samples_to_remove)

Hib_count_data <- Hib_count_data %>% purrr::set_names(c("S01B01", "S01B02", "S01B03", "S02B01", "S02B02", "S02B03", "S03B02", "S03B03", "S04B01", "S04B02", "S04B03", "S05B01", "S05B02", "S05B03", "S06B01", "S06B02", "S06B03"))

# Rename row names to old locus tag
# Rename Hib genes to old locus tag (available at KEGG)
# Sort row names alphabetically
# remove ribosomal RNA
Hib_count_data <- Hib_count_data %>%
  rownames_to_column() %>%
  left_join(RefSeq_gff.genes, by = c("rowname" = "Name")) %>%
  mutate(Gene_names = if_else(is.na(OldLocusTag), rowname, OldLocusTag)) %>%
  arrange(Gene_names) %>%
  dplyr::filter(!str_detect(Gene_names, pattern = "HIB_r")) %>%
  column_to_rownames(var = "Gene_names")

Hib_count_data <- Hib_count_data %>% dplyr::select(-c(rowname, OldLocusTag, product))

