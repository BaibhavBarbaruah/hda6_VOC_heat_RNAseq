############################################################
## 01_build_dds_from_gene_counts_and_sample_table.R
## - Reads:
##     counts/gene_counts_matrix.txt
##     counts/sample_table.csv
## - Cleans sample names
## - Creates DESeqDataSet `dds`
## - Runs DESeq()
## - Saves `dds_from_gene_counts_and_sample_table.RData`
############################################################

library(DESeq2)

## 0. Set working directory to the project folder ----------------------------

# Adjust this if your folder name changes, but for now:
setwd("F:/RNA seq")

counts_file   <- "counts/gene_counts_matrix.txt"
metadata_file <- "counts/sample_table.csv"

## 1. Read data --------------------------------------------------------------

counts <- read.delim(
  counts_file,
  row.names   = 1,
  check.names = FALSE
)

sample_info <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE
)

# sample_info columns: "X", "sample", "genotype", "treatment"
sample_info <- sample_info[, c("sample", "genotype", "treatment")]

## 2. Create a 'condition' column (e.g. Col0_NS, Col0_T2H, hda6_T2H) --------

sample_info$condition <- paste(
  sample_info$genotype,
  sample_info$treatment,
  sep = "_"
)
sample_info$condition <- factor(sample_info$condition)

## 3. Helper to clean sample names on both sides -----------------------------

clean_names <- function(x) {
  x <- basename(x)                    # drop any path like "align/..."
  x <- sub("\\.bam$", "", x)          # drop .bam
  x <- sub("\\.fastq\\.gz$", "", x)   # drop .fastq.gz (just in case)
  x
}

counts_samples_clean <- clean_names(colnames(counts))
sample_samples_clean <- clean_names(sample_info$sample)

## 4. Match metadata rows to counts columns ---------------------------------

missing_in_metadata <- setdiff(counts_samples_clean, sample_samples_clean)
if (length(missing_in_metadata) > 0) {
  stop("These samples are in counts but not in sample_table (after cleaning):\n",
       paste(missing_in_metadata, collapse = ", "), "\n",
       "Check for typos in counts/gene_counts_matrix.txt or counts/sample_table.csv.")
}

idx <- match(counts_samples_clean, sample_samples_clean)
if (any(is.na(idx))) {
  stop("Some samples could not be matched (NA indices).\n",
       "Check sample names in counts and sample_table.csv.")
}

sample_info_aligned <- sample_info[idx, , drop = FALSE]

# Put the true column names of counts as rownames of colData
rownames(sample_info_aligned) <- colnames(counts)

## 5. Create DESeq2 object and run DESeq ------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = sample_info_aligned,
  design    = ~ condition
)

keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]

dds <- DESeq(dds)

## 6. Save the DESeqDataSet object ------------------------------------------

save(dds, file = "dds_from_gene_counts_and_sample_table.RData")

message("Done! Saved dds as 'dds_from_gene_counts_and_sample_table.RData' in F:/RNA seq")

