################################################################################
# SARS-CoV-2 Variant Analysis in R
# Author: [Your Name]
# Date: October 2024
# Description: Downstream analysis of variants from nf-core/viralrecon pipeline
################################################################################

# ==============================================================================
# 1. SETUP AND PACKAGE INSTALLATION
# ==============================================================================

# Install required packages (run once)
# Uncomment the lines below if packages are not installed
# install.packages("tidyverse")
# install.packages("ggpubr")
# install.packages("RColorBrewer")

# Load required libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

# Set working directory (adjust to your project path)
# setwd("~/sars-cov2-analysis")

# Verify working directory
getwd()

# Create output directories
dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

# ==============================================================================
# 2. DATA LOADING AND EXPLORATION
# ==============================================================================

# Load variant data from viralrecon output
var <- read.csv("results/viralrecon/variants_long_table/combined_variants.csv")

# Basic data exploration
cat("\n=== DATA OVERVIEW ===\n")
cat("Dataset dimensions:", dim(var), "\n")
cat("Number of variants:", nrow(var), "\n")
cat("Number of variables:", ncol(var), "\n\n")

# Display first few rows
cat("First 6 rows:\n")
print(head(var))

# Display structure
cat("\n=== DATA STRUCTURE ===\n")
str(var)

# Summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
summary(var)

# Check column names
cat("\n=== COLUMN NAMES ===\n")
print(colnames(var))

# Check data types
cat("\n=== DATA TYPES ===\n")
cat("SAMPLE:", class(var$SAMPLE), "\n")
cat("CHROM:", class(var$CHROM), "\n")
cat("DP:", class(var$DP), "\n")

# ==============================================================================
# 3. DATA TRANSFORMATION
# ==============================================================================

# Convert to tibble for better printing
var_tb <- as_tibble(var)
cat("\n=== CONVERTED TO TIBBLE ===\n")
print(var_tb)

# Add log2-transformed depth
var_tb_log <- var_tb %>% mutate(DP_log2 = log2(DP))
cat("\n=== ADDED LOG2 DEPTH ===\n")
print(select(var_tb_log, SAMPLE, REF, ALT, DP, DP_log2) %>% head())

# ==============================================================================
# 4. VARIANT COUNTING AND SUMMARIZATION
# ==============================================================================

# Count variants per sample
cat("\n=== VARIANTS PER SAMPLE ===\n")
variant_counts <- var_tb %>% count(SAMPLE, sort = TRUE)
print(variant_counts)

# Count variants per gene per sample
cat("\n=== VARIANTS PER GENE PER SAMPLE (top 10) ===\n")
gene_variants <- var_tb %>% count(SAMPLE, GENE, sort = TRUE)
print(head(gene_variants, 10))

# Count variant effects
cat("\n=== VARIANT EFFECTS ===\n")
effect_counts <- var_tb %>% count(EFFECT, sort = TRUE)
print(effect_counts)

# Count effects per sample
cat("\n=== EFFECTS PER SAMPLE ===\n")
effect_per_sample <- var_tb %>% count(EFFECT, SAMPLE, sort = TRUE)
print(effect_per_sample)

# ==============================================================================
# 5. DEPTH STATISTICS
# ==============================================================================

# Overall depth statistics
cat("\n=== DEPTH STATISTICS (ALL SAMPLES) ===\n")
cat("Maximum DP:", max(var_tb$DP), "\n")
cat("Minimum DP:", min(var_tb$DP), "\n")
cat("Mean DP:", round(mean(var_tb$DP), 2), "\n")
cat("Median DP:", median(var_tb$DP), "\n")

# Depth statistics per sample
cat("\n=== DEPTH STATISTICS PER SAMPLE ===\n")
depth_stats <- var_tb %>% 
  group_by(SAMPLE) %>% 
  summarize(
    n_variants = n(),
    max_depth = max(DP),
    min_depth = min(DP),
    mean_depth = round(mean(DP), 2),
    median_depth = median(DP)
  )
print(depth_stats)

# Export depth statistics
write.csv(depth_stats, "tables/depth_statistics.csv", row.names = FALSE)
cat("\nDepth statistics saved to: tables/depth_statistics.csv\n")

# ==============================================================================
# 6. DATA FILTERING EXAMPLES
# ==============================================================================

# Filter for specific sample
cat("\n=== FILTERING FOR SRR13500958 ===\n")
srr_variants <- filter(var_tb, SAMPLE == "SRR13500958") %>%
  select(CHROM, POS, REF, ALT, DP)
print(head(srr_variants, 10))

# High-confidence variants (DP >= 500)
cat("\n=== HIGH-CONFIDENCE VARIANTS (DP >= 500) ===\n")
high_conf <- var_tb %>%
  filter(SAMPLE == "SRR13500958" & DP >= 500) %>%
  select(CHROM, POS, REF, ALT, DP)
print(high_conf)

# Very high-confidence variants (DP >= 1000)
cat("\n=== VERY HIGH-CONFIDENCE VARIANTS (DP >= 1000) ===\n")
very_high_conf <- var_tb %>%
  filter(SAMPLE == "SRR13500958" & DP >= 1000) %>%
  select(CHROM, POS, REF, ALT, DP)
print(very_high_conf)

# Stop codon variants
cat("\n=== STOP CODON VARIANTS ===\n")
stop_variants <- filter(var_tb, EFFECT %in% c("stop_lost", "stop_gained")) %>%
  select(SAMPLE, CHROM, GENE, EFFECT, POS, REF, ALT)
print(stop_variants)

# ==============================================================================
# 7. VISUALIZATIONS
# ==============================================================================

cat("\n=== GENERATING VISUALIZATIONS ===\n")

# -------------------------------
# Plot 1: Depth Distribution (Boxplot)
# -------------------------------
cat("Creating Plot 1: Depth distribution boxplot...\n")

p1 <- ggplot(data = var_tb, aes(x = SAMPLE, y = DP, fill = SAMPLE)) + 
  geom_boxplot() + 
  ylim(0, 10000) + 
  scale_fill_brewer(palette = "RdYlBu") +
  labs(
    title = "Sequencing Depth Distribution across Samples",
    x = "Sample ID",
    y = "Depth (DP)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("figures/depth_distribution_boxplot.png", plot = p1, 
       width = 10, height = 6, dpi = 300)

# -------------------------------
# Plot 2: Depth Distribution (Points)
# -------------------------------
cat("Creating Plot 2: Depth distribution points...\n")

p2 <- ggplot(data = var_tb, aes(x = SAMPLE, y = DP, fill = SAMPLE)) + 
  geom_point(shape = 21, size = 3, alpha = 0.7) + 
  ylim(0, 10000) + 
  scale_fill_brewer(palette = "RdYlBu") +
  labs(
    title = "Variant Depth Distribution (Individual Points)",
    x = "Sample ID",
    y = "Depth (DP)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("figures/depth_distribution_points.png", plot = p2, 
       width = 10, height = 6, dpi = 300)

# -------------------------------
# Plot 3: Depth per Chromosome (Boxplot + Facet)
# -------------------------------
cat("Creating Plot 3: Depth per chromosome...\n")

p3 <- ggplot(data = var_tb, aes(x = CHROM, y = DP, fill = SAMPLE)) + 
  geom_boxplot() +
  ylim(0, 10000) + 
  scale_fill_brewer(palette = "RdYlBu") + 
  labs(
    title = "Depth Distribution per Chromosome",
    x = "Chromosome",
    y = "Depth (DP)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  facet_grid(. ~ SAMPLE)

ggsave("figures/depth_per_chromosome.png", plot = p3, 
       width = 14, height = 6, dpi = 300)

# -------------------------------
# Plot 4: Depth per Chromosome (Violin + Boxplot)
# -------------------------------
cat("Creating Plot 4: Depth per chromosome (violin)...\n")

p4 <- ggplot(data = var_tb, aes(x = CHROM, y = DP, fill = SAMPLE)) + 
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  ylim(0, 10000) + 
  scale_fill_brewer(palette = "RdYlBu") + 
  labs(
    title = "Depth Distribution per Chromosome (Violin Plot)",
    x = "Chromosome",
    y = "Depth (DP)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  facet_grid(. ~ SAMPLE)

ggsave("figures/depth_per_chromosome_violin.png", plot = p4, 
       width = 14, height = 6, dpi = 300)

# -------------------------------
# Plot 5: Variant Effects per Sample
# -------------------------------
cat("Creating Plot 5: Variant effects...\n")

p5 <- ggplot(data = var_tb, aes(y = EFFECT, fill = SAMPLE)) + 
  geom_bar() +
  scale_fill_brewer(palette = "RdBu") + 
  labs(
    title = "Variant Effects per Sample",
    x = "Count",
    y = "Effect Type"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave("figures/variant_effects.png", plot = p5, 
       width = 10, height = 8, dpi = 300)

# -------------------------------
# Plot 6: Depth vs Genome Position
# -------------------------------
cat("Creating Plot 6: Depth vs position...\n")

p6 <- ggplot(data = var_tb, aes(x = POS, y = DP, fill = SAMPLE)) + 
  geom_point(shape = 21, size = 3, alpha = 0.7) +
  scale_fill_brewer(palette = "RdBu") + 
  labs(
    title = "Sequencing Depth across Genome Position",
    x = "Genome Position (bp)",
    y = "Depth (DP)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave("figures/depth_vs_position.png", plot = p6, 
       width = 12, height = 6, dpi = 300)

# -------------------------------
# Plot 7: Total Depth vs Alternative Allele Depth
# -------------------------------
cat("Creating Plot 7: DP vs ALT_DP...\n")

p7 <- ggplot(data = var_tb, aes(x = ALT_DP, y = DP, fill = SAMPLE)) +
  geom_point(shape = 21, size = 3, alpha = 0.7) +
  scale_fill_brewer(palette = "RdBu") +
  labs(
    title = "Total Depth vs Alternative Allele Depth",
    x = "Alternative Allele Depth (ALT_DP)",
    y = "Total Depth (DP)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave("figures/dp_vs_altdp.png", plot = p7, 
       width = 10, height = 6, dpi = 300)

# -------------------------------
# Plot 8: Log-scale Depth Distribution
# -------------------------------
cat("Creating Plot 8: Log-scale depth...\n")

p8 <- ggplot(data = var_tb, aes(x = SAMPLE, y = DP, fill = SAMPLE)) + 
  geom_boxplot() + 
  scale_y_log10() +
  scale_fill_brewer(palette = "RdYlBu") +
  labs(
    title = "Depth Distribution (Log Scale)",
    x = "Sample ID",
    y = "Depth (DP, log10)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("figures/depth_distribution_log.png", plot = p8, 
       width = 10, height = 6, dpi = 300)

# ==============================================================================
# 8. EXPORT SUMMARY TABLES
# ==============================================================================

cat("\n=== EXPORTING SUMMARY TABLES ===\n")

# Effect summary per sample
effect_summary <- var_tb %>%
  count(SAMPLE, EFFECT) %>%
  pivot_wider(names_from = EFFECT, values_from = n, values_fill = 0)

write.csv(effect_summary, "tables/effect_summary.csv", row.names = FALSE)
cat("Effect summary saved to: tables/effect_summary.csv\n")

# Gene variant counts
gene_summary <- var_tb %>%
  count(GENE, EFFECT, sort = TRUE)

write.csv(gene_summary, "tables/gene_variant_summary.csv", row.names = FALSE)
cat("Gene variant summary saved to: tables/gene_variant_summary.csv\n")

# Stop codon variants
write.csv(stop_variants, "tables/stop_codon_variants.csv", row.names = FALSE)
cat("Stop codon variants saved to: tables/stop_codon_variants.csv\n")

# ==============================================================================
# 9. GENERATE FINAL REPORT
# ==============================================================================

cat("\n=== GENERATING ANALYSIS REPORT ===\n")

# Create comprehensive summary
report <- list(
  total_variants = nrow(var_tb),
  samples = unique(var_tb$SAMPLE),
  n_samples = length(unique(var_tb$SAMPLE)),
  depth_stats = depth_stats,
  variant_counts = variant_counts,
  effect_counts = effect_counts,
  n_stop_variants = nrow(stop_variants)
)

# Save report as RDS (R object)
saveRDS(report, "tables/analysis_report.rds")
cat("Analysis report saved to: tables/analysis_report.rds\n")

# Print final summary
cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Total variants analyzed:", report$total_variants, "\n")
cat("Number of samples:", report$n_samples, "\n")
cat("Stop codon variants found:", report$n_stop_variants, "\n")
cat("\nAll figures saved to: figures/\n")
cat("All tables saved to: tables/\n")
cat("\n=== END OF ANALYSIS ===\n")

################################################################################
# END OF SCRIPT
################################################################################