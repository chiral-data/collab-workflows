# Node 2: DESeq2 Analysis
# Purpose: Create DESeq2 dataset for statistical analysis of differentially expressed genes
# Inputs: counts.csv (from featureCounts)
# Outputs: merged.csv, dds.rds, res.rds

library(DESeq2)
library(tidyverse)

# =============================================================================
# CONFIGURATION FROM ENVIRONMENT VARIABLES
# =============================================================================

control_condition <- Sys.getenv("PARAM_CONTROL_CONDITION", "normal")
treatment_condition <- Sys.getenv("PARAM_TREATMENT_CONDITION", "treatment")
alpha <- as.numeric(Sys.getenv("PARAM_ALPHA", "0.05"))
gene_id_column <- Sys.getenv("PARAM_GENE_ID_COLUMN", "Geneid")
count_columns_param <- Sys.getenv("PARAM_COUNT_COLUMNS", "")
conditions_param <- Sys.getenv("PARAM_CONDITIONS", "")

cat("DESeq2 Analysis Parameters:\n")
cat("  Control condition:", control_condition, "\n")
cat("  Treatment condition:", treatment_condition, "\n")
cat("  Alpha (significance threshold):", alpha, "\n")

# =============================================================================
# LOAD DATA
# =============================================================================

cat("\n--- Loading count data ---\n")

# Find counts file
count_files <- list.files(pattern = "counts\\.csv$", full.names = TRUE)
if (length(count_files) == 0) {
  stop("No counts.csv file found")
}
count_file <- count_files[1]
cat("Loading:", count_file, "\n")

count.file <- read.csv(count_file, header = TRUE)
cat("Loaded table with", nrow(count.file), "rows and", ncol(count.file), "columns\n")
cat("Column names:", paste(colnames(count.file), collapse = ", "), "\n")

# =============================================================================
# EXTRACT COUNT MATRIX
# =============================================================================

cat("\n--- Extracting count matrix ---\n")

# Auto-detect count columns (numeric columns that could be counts)
# Skip first few annotation columns
all_cols <- colnames(count.file)

# Try to find numeric columns with count-like data
# Usually count columns are the last columns or columns with sample names
numeric_cols <- sapply(count.file, function(x) is.numeric(x) && all(x >= 0, na.rm = TRUE))
potential_count_cols <- which(numeric_cols)

# Filter out common annotation columns
annotation_patterns <- c("Length", "Start", "End", "Chr", "Strand", "^X$", "^X\\.", "row")
for (pattern in annotation_patterns) {
  exclude <- grepl(pattern, all_cols, ignore.case = TRUE)
  potential_count_cols <- potential_count_cols[!exclude[potential_count_cols]]
}

if (length(potential_count_cols) < 2) {
  stop("Could not auto-detect count columns. Please specify PARAM_COUNT_COLUMNS.")
}

# Use detected count columns
count_col_indices <- potential_count_cols
cat("Detected count columns:", paste(all_cols[count_col_indices], collapse = ", "), "\n")

# Extract count matrix
count_matrix <- count.file[, count_col_indices]

# Set row names from gene ID column
if (gene_id_column %in% colnames(count.file)) {
  rownames(count_matrix) <- count.file[[gene_id_column]]
} else {
  # Try common alternatives
  alt_cols <- c("Geneid", "GeneID", "gene_id", "Gene", "gene")
  for (col in alt_cols) {
    if (col %in% colnames(count.file)) {
      rownames(count_matrix) <- count.file[[col]]
      cat("Using", col, "as gene ID column\n")
      break
    }
  }
}

# Remove rows with NA values
count_matrix <- count_matrix[complete.cases(count_matrix), ]

# Convert to integer matrix (required by DESeq2)
count_matrix <- round(as.matrix(count_matrix))
storage.mode(count_matrix) <- "integer"

cat("Count matrix dimensions:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")

# =============================================================================
# PREPARE METADATA
# =============================================================================

cat("\n--- Preparing sample metadata ---\n")

n_samples <- ncol(count_matrix)

# Determine conditions for each sample
if (conditions_param != "" && conditions_param != "null") {
  conditions <- strsplit(conditions_param, ",")[[1]]
  conditions <- trimws(conditions)
} else {
  # Default: alternate between conditions or use sample names
  if (n_samples == 2) {
    conditions <- c(control_condition, treatment_condition)
  } else {
    # Try to infer from column names
    conditions <- rep(c(control_condition, treatment_condition), length.out = n_samples)
  }
}

if (length(conditions) != n_samples) {
  warning("Number of conditions does not match number of samples. Using alternating pattern.")
  conditions <- rep(c(control_condition, treatment_condition), length.out = n_samples)
}

col_data <- data.frame(
  row.names = colnames(count_matrix),
  condition = factor(conditions)
)

cat("Sample metadata:\n")
print(col_data)

# =============================================================================
# CREATE DESEQ2 DATASET AND RUN ANALYSIS
# =============================================================================

cat("\n--- Running DESeq2 ---\n")

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = col_data,
  design = ~ condition
)

dds <- DESeq(dds)

cat("DESeq2 analysis complete\n")

# =============================================================================
# EXTRACT RESULTS
# =============================================================================

cat("\n--- Extracting results ---\n")

# Get contrast based on available conditions
available_conditions <- levels(col_data$condition)
cat("Available conditions:", paste(available_conditions, collapse = ", "), "\n")

# Results with default alpha
res <- results(dds, contrast = c("condition", treatment_condition, control_condition))

# Results with specified alpha
res_alpha <- results(dds,
                     contrast = c("condition", treatment_condition, control_condition),
                     alpha = alpha)

# Print summaries
cat("\n--- Results Summary (default alpha) ---\n")
summary(res)

cat("\n--- Results Summary (alpha =", alpha, ") ---\n")
summary(res_alpha)

# =============================================================================
# PRINCIPAL COMPONENT ANALYSIS
# =============================================================================

cat("\n--- Performing PCA ---\n")

# r-log transform data
rld <- rlog(dds, blind = FALSE)

# Extract PCA data
pca_data <- plotPCA(rld, intgroup = c("condition"), returnData = TRUE)
percent_var <- round(attr(pca_data, "percentVar") * 100, digits = 1)

cat("PC1 variance:", percent_var[1], "%\n")
cat("PC2 variance:", percent_var[2], "%\n")

# Save PCA data for plotting node
write.csv(pca_data, file = "pca_data.csv", row.names = FALSE)
cat("Saved: pca_data.csv\n")

# =============================================================================
# MERGE AND SAVE RESULTS
# =============================================================================

cat("\n--- Saving results ---\n")

# Convert DESeq2 results to a data frame
res_df <- as.data.frame(res_alpha)

# Add gene names as a column
res_df$Geneid <- rownames(res_df)

# Merge DESeq2 results with full featureCounts file
merged <- merge(count.file, res_df, by = "Geneid", all.x = TRUE)

# Save outputs
write.csv(merged, file = "merged.csv", row.names = FALSE)
cat("Saved: merged.csv\n")

saveRDS(dds, file = "dds.rds")
cat("Saved: dds.rds\n")

saveRDS(res, file = "res.rds")
cat("Saved: res.rds\n")

# Also save the results with alpha threshold
saveRDS(res_alpha, file = "res_alpha.rds")
cat("Saved: res_alpha.rds\n")

cat("\nDESeq2 analysis completed successfully!\n")
