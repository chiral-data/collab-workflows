# Node 3: Visualization Plots
# Purpose: Generate MA plot and Volcano plot for differential expression results
# Inputs: dds.rds, res.rds, merged.csv
# Outputs: ma_plot.png, volcano_plot.png, pca_plot.png, merged_with_regulation.csv

library(DESeq2)
library(tidyverse)
library(ashr)
library(dplyr)
library(ggrepel)

# Try to load ggrastr, fall back to regular geom_point if not available
use_ggrastr <- requireNamespace("ggrastr", quietly = TRUE)
if (use_ggrastr) {
  library(ggrastr)
}

# =============================================================================
# CONFIGURATION FROM ENVIRONMENT VARIABLES
# =============================================================================

alpha <- as.numeric(Sys.getenv("PARAM_ALPHA", "0.05"))
plot_width <- as.integer(Sys.getenv("PARAM_PLOT_WIDTH", "10"))
plot_height <- as.integer(Sys.getenv("PARAM_PLOT_HEIGHT", "8"))
plot_dpi <- as.integer(Sys.getenv("PARAM_PLOT_DPI", "300"))
up_color <- Sys.getenv("PARAM_UP_COLOR", "firebrick")
down_color <- Sys.getenv("PARAM_DOWN_COLOR", "dodgerblue4")
none_color <- Sys.getenv("PARAM_NONE_COLOR", "honeydew4")
highlight_color <- Sys.getenv("PARAM_HIGHLIGHT_COLOR", "gold1")

cat("Plot Parameters:\n")
cat("  Alpha:", alpha, "\n")
cat("  Dimensions:", plot_width, "x", plot_height, "inches\n")
cat("  DPI:", plot_dpi, "\n")

# =============================================================================
# LOAD DATA
# =============================================================================

cat("\n--- Loading data ---\n")

# Load DESeq2 objects (check inputs/ folder first)
dds_path <- if (file.exists("inputs/dds.rds")) "inputs/dds.rds" else "dds.rds"
res_path <- if (file.exists("inputs/res.rds")) "inputs/res.rds" else "res.rds"
merged_path <- if (file.exists("inputs/merged.csv")) "inputs/merged.csv" else "merged.csv"

dds <- readRDS(dds_path)
res <- readRDS(res_path)
cat("Loaded:", dds_path, ",", res_path, "\n")

# Load merged data
merged <- read.csv(merged_path)
cat("Loaded:", merged_path, "(", nrow(merged), "rows)\n")

# =============================================================================
# PERFORM SHRINKAGE FOR MA PLOT
# =============================================================================

cat("\n--- Performing log fold change shrinkage ---\n")

# Perform shrinkage using ashr
shrunk_MA <- lfcShrink(dds = dds, res = res, coef = 2, type = "ashr")

# Convert to data frame and add gene IDs
shrunk_MA <- as.data.frame(shrunk_MA)

# Detect gene ID column name from merged data
gene_col <- intersect(c("GeneID", "Geneid", "gene_id", "Gene"), colnames(merged))[1]
shrunk_MA[[gene_col]] <- rownames(shrunk_MA)

# Keep only shrunken LFC and gene ID
shrunk_MA <- dplyr::select(shrunk_MA, log2FoldChange, all_of(gene_col))
colnames(shrunk_MA)[1] <- "log2FoldChange.shrink"

cat("Shrinkage complete\n")

# =============================================================================
# MERGE WITH ORIGINAL DATA AND FLAG REGULATION
# =============================================================================

cat("\n--- Processing data for plotting ---\n")

# Merge with original data
merged <- dplyr::inner_join(merged, shrunk_MA, by = gene_col)

# Flag significant genes
merged$regulated <- merged$padj < alpha
merged$regulated[is.na(merged$regulated)] <- FALSE

# Add regulation direction
merged <- mutate(merged, Regulation = case_when(
  regulated & log2FoldChange > 0 ~ "Up",
  regulated & log2FoldChange < 0 ~ "Down",
  TRUE ~ "None"
))

# Count regulated genes
n_up <- sum(merged$Regulation == "Up", na.rm = TRUE)
n_down <- sum(merged$Regulation == "Down", na.rm = TRUE)
n_none <- sum(merged$Regulation == "None", na.rm = TRUE)

cat("Regulation summary:\n")
cat("  Upregulated:", n_up, "\n")
cat("  Downregulated:", n_down, "\n")
cat("  Not significant:", n_none, "\n")

# Save merged data with regulation info
write.csv(merged, file = "merged_with_regulation.csv", row.names = FALSE)
cat("Saved: merged_with_regulation.csv\n")

# =============================================================================
# DETERMINE LABEL COLUMN
# =============================================================================

# Try to find a good column for gene labels
label_cols <- c("Gene", "gene", "Name", "gene_name", "Symbol", "symbol")
label_col <- NULL
for (col in label_cols) {
  if (col %in% colnames(merged)) {
    label_col <- col
    break
  }
}

if (is.null(label_col)) {
  # Use GeneID as label (detect actual column name)
  gene_col <- intersect(c("GeneID", "Geneid", "gene_id"), colnames(merged))[1]
  merged$Label <- merged[[gene_col]]
  label_col <- "Label"
} else {
  merged$Label <- merged[[label_col]]
}

cat("Using", label_col, "for gene labels\n")

# =============================================================================
# MA PLOT
# =============================================================================

cat("\n--- Generating MA plot ---\n")

# Create base plot
ma_plot <- ggplot(merged, aes(x = log10(baseMean),
                              y = log2FoldChange.shrink,
                              col = Regulation,
                              label = Label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("Up" = up_color, "Down" = down_color, "None" = none_color)) +
  theme_classic() +
  labs(
    title = "MA Plot",
    subtitle = paste0("Upregulated: ", n_up, " | Downregulated: ", n_down, " (padj < ", alpha, ")"),
    x = "log10 Base Mean",
    y = "log2 Fold Change (shrunk)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "right"
  )

# Add points (use ggrastr if available for performance)
if (use_ggrastr) {
  ma_plot <- ma_plot + ggrastr::rasterise(geom_point(alpha = 0.6), dpi = plot_dpi)
} else {
  ma_plot <- ma_plot + geom_point(alpha = 0.6)
}

# Add labels for genes with names
if (any(!is.na(merged$Label) & merged$Label != "")) {
  labeled_genes <- subset(merged, !is.na(Label) & Label != "" & regulated)
  if (nrow(labeled_genes) > 0) {
    # Limit to top genes by significance
    labeled_genes <- labeled_genes[order(labeled_genes$padj), ]
    labeled_genes <- head(labeled_genes, 20)

    ma_plot <- ma_plot +
      geom_point(data = labeled_genes, col = highlight_color, size = 2) +
      geom_label_repel(data = labeled_genes, show.legend = FALSE,
                       max.overlaps = 15, size = 3)
  }
}

# Save MA plot
tryCatch({
  ggsave("ma_plot.png", plot = ma_plot, width = plot_width, height = plot_height, dpi = plot_dpi)
  cat("Saved: ma_plot.png\n")
}, error = function(e) {
  cat("Error saving MA plot:", conditionMessage(e), "\n")
})

# =============================================================================
# VOLCANO PLOT
# =============================================================================

cat("\n--- Generating Volcano plot ---\n")

# Create base plot
volcano_plot <- ggplot(merged, aes(x = log2FoldChange.shrink,
                                    y = -log10(padj),
                                    col = Regulation,
                                    label = Label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("Up" = up_color, "Down" = down_color, "None" = none_color)) +
  theme_classic() +
  labs(
    title = "Volcano Plot",
    subtitle = paste0("Upregulated: ", n_up, " | Downregulated: ", n_down, " (padj < ", alpha, ")"),
    x = "log2 Fold Change (shrunk)",
    y = "-log10 Adjusted p-value"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "right"
  )

# Add points
if (use_ggrastr) {
  volcano_plot <- volcano_plot + ggrastr::rasterise(geom_point(alpha = 0.6), dpi = plot_dpi)
} else {
  volcano_plot <- volcano_plot + geom_point(alpha = 0.6)
}

# Add labels for genes with names
if (any(!is.na(merged$Label) & merged$Label != "")) {
  labeled_genes <- subset(merged, !is.na(Label) & Label != "" & regulated)
  if (nrow(labeled_genes) > 0) {
    labeled_genes <- labeled_genes[order(labeled_genes$padj), ]
    labeled_genes <- head(labeled_genes, 20)

    volcano_plot <- volcano_plot +
      geom_point(data = labeled_genes, col = highlight_color, size = 2) +
      geom_label_repel(data = labeled_genes, show.legend = FALSE,
                       max.overlaps = 15, size = 3)
  }
}

# Save Volcano plot
tryCatch({
  ggsave("volcano_plot.png", plot = volcano_plot, width = plot_width, height = plot_height, dpi = plot_dpi)
  cat("Saved: volcano_plot.png\n")
}, error = function(e) {
  cat("Error saving Volcano plot:", conditionMessage(e), "\n")
})

# =============================================================================
# PCA PLOT
# =============================================================================

cat("\n--- Generating PCA plot ---\n")

# Check if PCA data exists (check inputs/ folder first)
pca_path <- if (file.exists("inputs/pca_data.csv")) "inputs/pca_data.csv" else "pca_data.csv"
if (file.exists(pca_path)) {
  pca_data <- read.csv(pca_path)

  percent_var <- c(unique(pca_data$PC1_var), unique(pca_data$PC2_var))
  
  # Create PCA plot
  pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, colour = condition)) +
    geom_point(size = 4) +
    theme_classic() +
    labs(
      title = "PCA Plot",
      x = paste0("PC1: ", percent_var[1], "% variance"),
      y = paste0("PC2: ", percent_var[2], "% variance")
    ) +
    scale_color_manual(values = c('#E69F00', '#56B4E9')) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      aspect.ratio = 1,
      legend.position = "right"
    )

  ggsave("pca_plot.png", plot = pca_plot, width = plot_height, height = plot_height, dpi = plot_dpi)
  cat("Saved: pca_plot.png\n")
} else {
  cat("PCA data not found, skipping PCA plot\n")
}

cat("\nPlot generation completed successfully!\n")
