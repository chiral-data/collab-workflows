# Node 1: featureCounts
# Purpose: Quantify mapped reads per gene using Rsubread's featureCounts
# Inputs: BAM files (aligned reads), GTF file (genome annotation)
# Outputs: counts.csv, counts_annotation.csv, counts_log.txt

library(Rsubread)

# =============================================================================
# CONFIGURATION FROM ENVIRONMENT VARIABLES
# =============================================================================

# Get parameters from environment
threads <- as.integer(Sys.getenv("PARAM_THREADS", "8"))
strandedness <- as.integer(Sys.getenv("PARAM_STRANDEDNESS", "1"))
feature_type <- Sys.getenv("PARAM_FEATURE_TYPE", "gene")
attribute_type <- Sys.getenv("PARAM_ATTRIBUTE_TYPE", "locus_tag")
long_read <- tolower(Sys.getenv("PARAM_LONG_READ", "false")) == "true"
allow_multi_overlap <- tolower(Sys.getenv("PARAM_ALLOW_MULTI_OVERLAP", "true")) == "true"
use_fraction <- tolower(Sys.getenv("PARAM_USE_FRACTION", "true")) == "true"
extra_attributes <- Sys.getenv("PARAM_EXTRA_ATTRIBUTES",
                               "ID,gene_biotype,Dbxref,locus_tag,Name,ECOCYC,gene_synonym,gene,orig_transcript_id")

# =============================================================================
# FIND INPUT FILES
# =============================================================================

# Find BAM files
bam_files <- list.files(pattern = "\\.bam$", full.names = TRUE)
if (length(bam_files) == 0) {
  stop("No BAM files found in the current directory")
}
cat("Found BAM files:\n")
print(bam_files)

# Find GTF file
gtf_files <- list.files(pattern = "\\.gtf$", full.names = TRUE)
if (length(gtf_files) == 0) {
  stop("No GTF file found in the current directory")
}
gtf_file <- gtf_files[1]
cat("Using GTF file:", gtf_file, "\n")

# =============================================================================
# RUN FEATURECOUNTS
# =============================================================================

cat("\n--- Running featureCounts ---\n")
cat("Parameters:\n")
cat("  Threads:", threads, "\n")
cat("  Strandedness:", strandedness, "(0=unstranded, 1=stranded, 2=reverse)\n")
cat("  Feature type:", feature_type, "\n")
cat("  Attribute type:", attribute_type, "\n")
cat("  Long read mode:", long_read, "\n")
cat("  Allow multi-overlap:", allow_multi_overlap, "\n")
cat("  Fractional counting:", use_fraction, "\n")

# Parse extra attributes
extra_cols <- strsplit(extra_attributes, ",")[[1]]
extra_cols <- trimws(extra_cols)

# Capture output for log
log_con <- file("counts_log.txt", open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

fc <- featureCounts(
  files = bam_files,
  annot.ext = gtf_file,
  isGTFAnnotationFile = TRUE,

  # Threads
  nthreads = threads,

  # Strandedness
  strandSpecific = strandedness,

  # Feature & attribute
  GTF.featureType = feature_type,
  GTF.attrType = attribute_type,

  # Feature-level counting (not meta-features)
  useMetaFeatures = FALSE,

  # Overlapping + fractional counts
  allowMultiOverlap = allow_multi_overlap,
  fraction = use_fraction,

  # Primary alignments only
  primaryOnly = TRUE,

  # Long read mode if specified
  isLongRead = long_read
)

sink(type = "message")
sink(type = "output")
close(log_con)

# =============================================================================
# SAVE OUTPUTS
# =============================================================================

cat("\n--- Saving outputs ---\n")

# Save counts matrix
write.csv(fc$counts, file = "counts_table.csv")
cat("Saved: counts_table.csv\n")

# Save annotation
write.csv(fc$annotation, file = "counts_annotation.csv")
cat("Saved: counts_annotation.csv\n")

# Create combined output (annotation + counts)
combined_table <- cbind(fc$annotation, fc$counts)
write.csv(combined_table, file = "counts.csv", row.names = FALSE)
cat("Saved: counts.csv (combined annotation + counts)\n")

# Print summary statistics
cat("\n--- Summary ---\n")
cat("Total features counted:", nrow(fc$counts), "\n")
cat("Number of samples:", ncol(fc$counts), "\n")
cat("Sample names:", paste(colnames(fc$counts), collapse = ", "), "\n")

cat("\nfeatureCounts completed successfully!\n")
