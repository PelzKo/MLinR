#!/usr/bin/env Rscript
# ==============================================================================
# Prepare_COVID_Data.R
# Data preparation script for the COVID-19 ML project
#
# This script downloads and processes the GSE157103 dataset for use in the
# "Machine Learning with R in the Life Sciences" course (Day 2 project).
#
# Source: Overmyer et al. (2021) Cell Systems 12(1):23-40
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157103
#
# The dataset contains bulk RNA-seq of peripheral blood leukocytes from
# 126 hospitalized patients:
#   - 100 COVID-19 positive (50 ICU, 50 non-ICU)
#   -  26 non-COVID-19      (16 ICU, 10 non-ICU)
#
# Output files:
#   - covid_expression.csv: log2(TPM+1) expression of top 5000 variable genes
#   - covid_metadata.csv:   clinical metadata for all 126 samples
# ==============================================================================

library(tidyverse)
#library(GEOquery)

# ==============================================================================
# 1. Download data from GEO
# ==============================================================================

# The TPM and expected count matrices are available as supplementary files.
# Download them manually from GEO or use the URLs below.

# Option A: Download manually from GEO
# Go to: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157103
# Download: GSE157103_genes.tpm.tsv.gz
#           GSE157103_series_matrix.txt.gz

# Option B: Download programmatically (may require GEOquery)
# gse <- getGEO("GSE157103", GSEMatrix = TRUE)
# series_matrix <- gse[[1]]

# For this script, we assume the files are already downloaded and in
# the working directory.

cat("=== COVID-19 Data Preparation ===\n\n")

# ==============================================================================
# 2. Parse metadata from series matrix
# ==============================================================================

cat("Step 1: Parsing metadata...\n")

# Read the series matrix file
# The series matrix is a special GEO format with metadata in comment lines
series_lines <- read_lines("GSE157103_series_matrix.txt.gz")

# Extract sample characteristics
# These lines start with !Sample_characteristics_ch1
char_lines <- series_lines[str_detect(series_lines, "^!Sample_characteristics_ch1")]

# Extract sample IDs and GEO accessions
geo_ids <- series_lines[str_detect(series_lines, "^!Sample_geo_accession")] |>
  str_remove("^!Sample_geo_accession\t") |>
  str_split("\t") |>
  pluck(1) |>
  str_remove_all('"')

sample_titles <- series_lines[str_detect(series_lines, "^!Sample_title")] |>
  str_remove("^!Sample_title\t") |>
  str_split("\t") |>
  pluck(1) |>
  str_remove_all('"')

# Parse each characteristic line into a tidy format
parse_characteristics <- function(char_lines) {
  result <- list()

  for (line in char_lines) {
    parts <- str_split(line, "\t")[[1]]
    parts <- str_remove_all(parts[-1], '"')  # Remove header and quotes

    # Each entry is "key: value"
    key <- str_extract(parts[1], "^[^:]+") |> str_trim()
    values <- str_remove(parts, "^[^:]+:\\s*")

    # Clean key name for column use
    clean_key <- key |>
      str_to_lower() |>
      str_replace_all("[^a-z0-9]+", "_") |>
      str_remove("_$")

    result[[clean_key]] <- values
  }

  as_tibble(result)
}

meta_raw <- parse_characteristics(char_lines)

# Extract sample IDs from titles (e.g., "COVID_01_39y_male_NonICU" -> "C1")
# Sample naming: C1-C100 for COVID, NC1-NC26 for non-COVID
# This matches the column names in the expression matrix

# Read the TPM file header to get exact sample IDs
tpm_header <- read_lines("GSE157103_genes.tpm.tsv.gz", n_max = 1) |>
  str_split("\t") |>
  pluck(1)
sample_ids <- tpm_header[-1]  # Remove gene name column

# Build metadata table
meta <- tibble(
  sample_id = sample_ids,
  geo_accession = geo_ids,
  sample_title = sample_titles
) |>
  bind_cols(meta_raw)

# Standardize column names and missing values
meta <- meta |>
  rename_with(~ str_replace_all(.x, "[^a-z0-9_]", "_")) |>
  mutate(across(everything(), ~ if_else(.x %in% c("", "unknown", "Unknown", "N/A"), NA_character_, .x)))

cat("  Metadata parsed:", nrow(meta), "samples,", ncol(meta), "variables\n")
cat("  Disease states:", paste(table(meta$disease_state), collapse = ", "), "\n")

# ==============================================================================
# 3. Process TPM expression data
# ==============================================================================

cat("\nStep 2: Processing expression data...\n")

# Read the full TPM matrix
tpm <- read_tsv("GSE157103_genes.tpm.tsv.gz", show_col_types = FALSE)
gene_names <- tpm[[1]]
tpm_matrix <- tpm[, -1] |> as.matrix()
rownames(tpm_matrix) <- gene_names

cat("  Raw matrix:", nrow(tpm_matrix), "genes x", ncol(tpm_matrix), "samples\n")

# ==============================================================================
# 4. Filter genes
# ==============================================================================

cat("\nStep 3: Filtering genes...\n")

# Filter 1: Require expression in at least 10% of samples (TPM > 1)
n_samples <- ncol(tpm_matrix)
expressed <- rowSums(tpm_matrix > 1) >= (n_samples * 0.10)

cat("  Genes with >=10% samples TPM>1:", sum(expressed), "\n")

# Filter 2: Require non-zero variance on log2 scale
log_tpm <- log2(tpm_matrix[expressed, ] + 1)
gene_vars <- apply(log_tpm, 1, var)
nonzero_var <- gene_vars > 0

cat("  Genes with non-zero variance:", sum(nonzero_var), "\n")

# ==============================================================================
# 5. Select top 5000 most variable genes
# ==============================================================================

cat("\nStep 4: Selecting top 5000 most variable genes...\n")

gene_vars_filtered <- gene_vars[nonzero_var]
top_5000_names <- names(sort(gene_vars_filtered, decreasing = TRUE)[1:5000])

# Log2(TPM+1) transform the selected genes
expr_final <- log2(tpm_matrix[top_5000_names, ] + 1)

# Transpose to samples-as-rows format
expr_df <- as_tibble(t(expr_final)) |>
  mutate(sample_id = sample_ids, .before = 1)

cat("  Final expression matrix:", nrow(expr_df), "samples x",
    ncol(expr_df) - 1, "genes\n")

# Check key COVID genes
key_genes <- c("IFI27", "ISG15", "MX1", "OAS1", "S100A8", "S100A12",
               "IFITM3", "IL1B", "CCL2", "DEFA3")
cat("  Key COVID genes in top 5000:\n")
for (g in key_genes) {
  if (g %in% top_5000_names) {
    rank <- which(names(sort(gene_vars_filtered, decreasing = TRUE)) == g)
    cat("    ", g, ": YES (variance rank", rank, ")\n")
  } else {
    cat("    ", g, ": no\n")
  }
}

# ==============================================================================
# 6. Write output files
# ==============================================================================

cat("\nStep 5: Writing output files...\n")

# Round expression values for manageable file size
expr_df <- expr_df |>
  mutate(across(-sample_id, ~ round(.x, 4)))

write_csv(expr_df, "covid_expression.csv")
write_csv(meta, "covid_metadata.csv")

cat("  covid_expression.csv:", file.size("covid_expression.csv") / 1e6, "MB\n")
cat("  covid_metadata.csv:", file.size("covid_metadata.csv") / 1e3, "KB\n")

# ==============================================================================
# 7. Summary
# ==============================================================================

cat("\n=== Summary ===\n")
cat("Dataset: GSE157103 (Overmyer et al. 2021)\n")
cat("Samples:", nrow(meta), "\n")
cat("  COVID-19:", sum(meta$disease_state == "COVID-19"), "\n")
cat("    ICU:", sum(meta$disease_state == "COVID-19" & meta$icu == "yes", na.rm = TRUE), "\n")
cat("    Non-ICU:", sum(meta$disease_state == "COVID-19" & meta$icu == "no", na.rm = TRUE), "\n")
cat("  Non-COVID:", sum(meta$disease_state == "non-COVID-19"), "\n")
cat("Genes: 5000 (top variable, log2(TPM+1) transformed)\n")
cat("Files: covid_expression.csv, covid_metadata.csv\n")
cat("\nDone!\n")
