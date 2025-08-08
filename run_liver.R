# Load custom functions from external file
source("functions.R")

# Set directory path for liver data
dir <- "human_adult_liver"

# Read miRNA-gene interaction data for liver
mirna_gene_liver <- read_csv(file.path(dir, "mirna_gene_organ_valid.csv"))

# Read gene-miRNA interaction data for liver and preprocess
gene_mirna_liver <- read_csv(file.path(dir, "gene_mirna_organ_valid.csv")) %>%
  
  # Remove genes not found in SCE
  filter(VG != "SCE NotFound") %>%
  
  # Select relevant columns
  dplyr::select(Gene.ID, target.number, miRNA)

# Read gene statistics and join with miRNA target information
gene_stats <- read_csv(file.path(dir, "gene_stats.csv")) %>%
  
  # Filter genes with UTR sequence available
  filter(is.na(UTR_sequence) == FALSE) %>%
  
  # Join with miRNA data
  left_join(gene_mirna_liver, by = c("ensembl_gene_id" = "Gene.ID"))

# Create subset of genes that are not miRNA targets
gene_stats_nt <- gene_stats %>%
  filter(is_target == FALSE)



# Create directory for saving statistical analysis plots if it doesn't exist
sa_path <- file.path(dir, "Statistical Analysis plots")
if (!dir.exists(sa_path)) {
  dir.create(sa_path, recursive = TRUE)
}

# Create Q-Q plot of gene expression noise (Delta)
p1 <- ggplot(gene_stats, aes(sample = log(Delta))) +
  stat_qq() +
  stat_qq_line(color = "red") +
  labs(title = "Q-Q Plot (Liver)", 
       x = "Theoretical Quantiles", 
       y = "Sample Quantiles") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )
p1_path <- file.path(sa_path, "QQ Plot.png")
ggsave(p1_path, p1, width = 12, height = 6, dpi = 600)



# Function to perform Wilcoxon test and create comparison plot for different gene groups
wilcox_test_vg <- function(vg) {
  
  # Validate input parameter
  if (!vg %in% c("HVG", "LVG", "Neither", "All")) {
    stop("vg must be one of: 'HVG', 'LVG', 'Neither', or 'All'")
  }
  
  # Filter data based on gene group
  if (vg == "HVG") {
    gene_stats_filtered <- gene_stats %>% filter(VG == "HVG")
  } else if (vg == "LVG") {
    gene_stats_filtered <- gene_stats %>% filter(VG == "LVG")
  } else if (vg == "Neither") {
    gene_stats_filtered <- gene_stats %>% filter(VG == "Not HVG or LVG")
  } else { # "All"
    gene_stats_filtered <- gene_stats
  }
  
  # Perform Wilcoxon test
  wilcox.test(log(Delta) ~ miRNA_target, data = gene_stats_filtered)
  mean_comparison <- wilcox.test(
    log(Delta) ~ miRNA_target, 
    data = gene_stats_filtered
  )
  
  # Format p-value and W statistic for display
  p_label <- ifelse(
    mean_comparison$p.value < 0.001,
    "p < 0.001", 
    sprintf("p = %.3f", mean_comparison$p.value)
  )
  w_label <- sprintf("w = %.0f", mean_comparison$statistic)
  
  # Create violin-boxplot comparing expression noise between target and non-target genes
  violin_box <- ggplot(
    gene_stats_filtered,
    aes(x = miRNA_target, y = log(Delta), fill = miRNA_target)
  ) +
    geom_violin(alpha = 0.5, width = 0.7) +
    geom_boxplot(
      width = 0.2,
      outlier.shape = NA,
      alpha = 0.7
    ) +
    scale_fill_manual(values = c("skyblue", "salmon")) +
    labs(
      x = NULL, y = "log-scaled Gene Expression Noise",
      title = paste0("Liver Gene Expression Noise Comparison (", vg, ")"),
      subtitle = paste("Wilcoxon test:", w_label, ",", p_label)
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "none"
    )
  violin_box
}

# Create and save plots for different gene groups
p2 <- wilcox_test_vg("HVG")
p2_path <- file.path(sa_path, "HVG.png")
ggsave(p2_path, p2, width = 8, height = 6, dpi = 600)

p3 <- wilcox_test_vg("LVG")
p3_path <- file.path(sa_path, "LVG.png")
ggsave(p3_path, p3, width = 8, height = 6, dpi = 600)

p4 <- wilcox_test_vg("All")
p4_path <- file.path(sa_path, "All.png")
ggsave(p4_path, p4, width = 8, height = 6, dpi = 600)



# Plot distribution of miRNA target numbers per gene
p5 <- ggplot(gene_stats,
             aes(x = target.number)) +
  geom_histogram(fill = "darkgreen", alpha = 0.8, color = "white") +
  labs(title = "Distribution of miRNA Target Numbers in Liver Genes",
       x = "Number of miRNA target",
       y = "Count") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )
p5_path <- file.path(sa_path, "miRNA distribution.png")
ggsave(p5_path, p5, width = 8, height = 6, dpi = 600)



# Compare expression noise between genes with highest and lowest miRNA target numbers
gene_stats_mirna <- gene_stats %>% 
  filter(is_target == TRUE)
quantiles <- quantile(gene_stats_mirna$target.number, probs = c(0.1, 0.9), na.rm = TRUE)
gene_stats_mirna <- gene_stats_mirna %>%
  mutate(
    target_group = case_when(
      target.number >= quantiles[2] ~ "Top10%",
      target.number <= quantiles[1] ~ "Bottom10%",
      TRUE ~ "Middle80" 
    )
  )
gene_stats_comparison <- gene_stats_mirna %>%
  filter(target_group %in% c("Top10%", "Bottom10%"))

# Perform Wilcoxon test between top and bottom groups
wilcox_result <- wilcox.test(
  log(Delta) ~ target_group,
  data = gene_stats_comparison
)

# Format p-value and W statistic for display
p_value <- wilcox_result$p.value
w_statistic <- wilcox_result$statistic
p_label <- ifelse(
  p_value < 0.001,
  "p < 0.001",
  sprintf("p = %.3f", p_value)
)
w_label <- sprintf("W = %.0f", w_statistic)

# Create comparison plot
p6 <- ggplot(
  gene_stats_comparison,
  aes(x = target_group, y = log(Delta), fill = target_group)
) +
  geom_violin(alpha = 0.5, width = 0.7) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.7) +
  scale_fill_manual(values = c("skyblue", "salmon")) +
  labs(
    x = "Target Number Group",
    y = "log-scaled Gene Expression Noise",
    title = "Expression Noise Comparison by miRNA Target Number (Liver)",
    subtitle = paste("Wilcoxon test:", w_label, ",", p_label)
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "none"
  )
p6_path <- file.path(sa_path, "Target Number Comparison.png")
ggsave(p6_path, p6, width = 8, height = 6, dpi = 600)



# Identify most highly expressed miRNA isoforms
mirna_high_expr <- mirna_gene_liver %>%
  arrange(desc(exp.mean)) %>%
  
  # Remove isoform suffixes
  mutate(miR_prefix = str_replace(miRNA, "-[a-z]?[35]p$", "")) %>%
  
  # Group by miRNA family
  group_by(miR_prefix) %>%
  
  # Take highest expressed isoform
  slice_max(order_by = exp.mean, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  pull(miRNA)

# Filter miRNA-gene data to only include most highly expressed isoforms
mirna_gene_filtered_liver <- mirna_gene_liver %>%
  filter(miRNA %in% mirna_high_expr) %>%
  arrange(desc(exp.mean))

# Get top 5 most highly expressed miRNAs and save to file
mirna_gene_liver_5 <- mirna_gene_filtered_liver[1:5,]
top_mir_path <-  file.path(sa_path, "top_five_mirna.csv")
write_csv(mirna_gene_liver_5, top_mir_path)



# Function to check proportion of target genes for a specific miRNA
check_proportion <- function(mir) {
  
  # Expand miRNA-gene data (one row per gene-miRNA pair)
  many_many <- mirna_gene_filtered_liver %>%
    dplyr::select(-gene.number) %>%
    separate_rows(Gene.ID, sep = "/")
  
  # Get target genes for specified miRNA
  mir <- many_many %>% filter(miRNA == mir) %>% pull(Gene.ID)
  
  # Get statistics for target genes (keeping longest UTR per gene)
  gene_stats_mir <- gene_stats[gene_stats$ensembl_gene_id %in% mir, ] %>%
    group_by(ensembl_gene_id) %>%
    slice_max(utr_length, n = 1) %>%
    ungroup()
  
  # Summarize statistics by gene group and target status
  stats_summary_mir <- gene_stats_mir %>%
    group_by(ensembl_gene_id) %>%
    slice_max(utr_length, n = 1) %>%
    ungroup() %>%
    group_by(VG, is_target) %>%
    summarise(
      mean_utr_length = mean(utr_length, na.rm = TRUE),
      count = n(),
      UTR_missing = sum(is.na(utr_length))
    ) %>%
    ungroup() %>%
    arrange(VG, is_target) %>%
    mutate(
      is_target = ifelse(is_target, "Target", "Non-target")
    )
  return(list(gene_stats_mir = gene_stats_mir,
              stats_summary_mir = stats_summary_mir))
}



# Convert data and add logarithmic UTR lengths
gene_stats_ecdf <- gene_stats %>% 
  filter(VG != "Not HVG or LVG") %>%
  mutate(log_utr = log(utr_length + 1))  # +1 to avoid taking logarithms to 0

wilcox_result <- wilcox.test(log_utr ~ VG, data = gene_stats_ecdf)
p_label <- ifelse(wilcox_result$p.value < 0.001, 
                  "p < 0.001", 
                  sprintf("p = %.3f", wilcox_result$p.value))
w_label <- sprintf("W = %.0f", wilcox_result$statistic)

p7 <- ggplot(gene_stats_ecdf, 
             aes(x = log_utr, color = VG)) +
  stat_ecdf(linewidth = 1, alpha = 0.8) + 
  labs(x = "Log-scaled UTR Length", 
       y = "Cumulative Probability",
       title = "eCDF of UTR Length by Gene Type (Liver)",
       subtitle = paste("Wilcoxon test:", w_label, ",", p_label)) +
  scale_color_manual(values = c("salmon", "skyblue"), 
                     name = "") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "right"
  )
p7_path <- file.path(sa_path, "UTR Length comparison.png")
ggsave(p7_path, p7, width = 8, height = 6, dpi = 600)



# Function to analyze individual miRNA effects on gene expression noise
individual_mir_check <- function(mir, vg = NULL) {
  
  # Get target gene statistics for specified miRNA
  check <- check_proportion(mir)
  gene_stats_mir <- check$gene_stats_mir
  stats_summary_mir <- check$stats_summary_mir
  
  # Filter data by gene group if specified
  if (is.null(vg)) {
    gene_stats_nt_vg <- gene_stats_nt
    gene_stats_mir_vg <- gene_stats_mir
    plot_title_suffix <- "(All)"
  } else {
    gene_stats_nt_vg <- gene_stats_nt %>% filter(VG == vg)
    gene_stats_mir_vg <- gene_stats_mir %>% filter(VG == vg)
    plot_title_suffix <- paste("(", vg, ")")
  }
  
  # Get expression noise values for comparison
  mir_vg_delta <- gene_stats_mir_vg$Delta
  non_target_vg_delta <- gene_stats_nt_vg$Delta
  
  # Perform Wilcoxon test between target and non-target genes
  wilcox_result <- wilcox.test(
    x = mir_vg_delta,
    y = non_target_vg_delta,
    alternative = "two.sided",
    paired = FALSE,
    exact = FALSE,
    conf.int = TRUE
  )
  
  # Format p-value and W statistic for display
  p_label <- ifelse(
    wilcox_result$p.value < 0.001,
    "p < 0.001",
    sprintf("p = %.3f", wilcox_result$p.value)
  )
  w_label <- sprintf("W = %.0f", wilcox_result$statistic)
  
  # Prepare data for plotting
  combined_data <- data.frame(
    Delta = c(non_target_vg_delta, mir_vg_delta),
    Group = rep(c("Non-targets", paste0(mir, " Targets")), 
                times = c(length(non_target_vg_delta), length(mir_vg_delta)))
  )
  
  # Create comparison plot
  violin_box <- ggplot(
    combined_data,
    aes(x = Group, y = log(Delta), fill = Group)
  ) +
    geom_violin(alpha = 0.5, width = 0.7) +
    geom_boxplot(
      width = 0.2,
      outlier.shape = NA,
      alpha = 0.7
    ) +
    scale_fill_manual(values = c("salmon", "skyblue")) +
    labs(
      x = NULL, 
      y = "log-scaled Gene Expression Noise",
      title = paste("Liver Gene Expression Noise Comparison", plot_title_suffix),
      subtitle = paste("Wilcoxon test:", w_label, ",", p_label)
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "none"
    )
  
  return(list(gene_stats_mir = gene_stats_mir, violin_box = violin_box))
}

# Analyze and plot effect of miR-122-5p
p8 <- individual_mir_check("miR-122-5p")$violinbox
p8_path <- file.path(sa_path, "mir_122_5p_check.png")
ggsave(p8_path, p8, width = 8, height = 6, dpi = 600)



# Function to perform GO enrichment analysis for miRNA target genes
run_go_enrichment <- function(mir,
                              vg = NULL,
                              organism = "org.Hs.eg.db",
                              ont = "BP",
                              pvalue_cutoff = 0.05,
                              qvalue_cutoff = 0.1,
                              show_category = 10) {
  
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(enrichplot)
    library(organism, character.only = TRUE)
  })
  
  # Get target gene statistics for specified miRNA
  gene_stats_mir <- individual_mir_check(mir, vg)$gene_stats_mir
  miR_targets <- gene_stats_mir %>% pull(ensembl_gene_id)
  
  # Perform GO enrichment analysis
  ego <- enrichGO(
    gene = miR_targets,
    OrgDb = eval(parse(text = organism)),
    keyType = "ENSEMBL",
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = pvalue_cutoff,
    qvalueCutoff = qvalue_cutoff,
    readable = TRUE
  )
  
  plots <- list()
  if (nrow(ego) > 0) {
    
    # Create bar plot of enriched GO terms
    plots$barplot <- barplot(
      ego,
      showCategory = show_category
    ) +
      labs(
        x = "Gene Count",
        title = paste("Liver ", mir, " GO Enrichment")
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
      )
    
    # Create network plot if enough terms are enriched
    if (nrow(ego) >= 3) {
      plots$cnetplot <- cnetplot(
        ego, 
        categorySize = "pvalue",
        showCategory = min(show_category, nrow(ego)),
        max.overlaps = 50)
    }
  } else {
    warning("No terms enriched under current cutoff!")
  }
  
  return(list(
    enrich_result = ego,
    plots = plots
  ))
}

# Perform GO enrichment for miR-122-5p and create plots
result1 <- run_go_enrichment("miR-122-5p", pvalue_cutoff = 0.1)

p9 <- result1$plots$barplot
p9_path <- file.path(sa_path, "GO barplot miR122.png")
ggsave(p9_path, p9, width = 8, height = 6, dpi = 600)

p10 <- result1$plots$cnetplot
p10_path <- file.path(sa_path, "GO cnetplot miR122.png")
ggsave(p10_path, p10, width = 8, height = 6, dpi = 600)

# Perform GO enrichment for miR-26a-5p and create plots
result2 <- run_go_enrichment("miR-26a-5p", pvalue_cutoff = 0.1)

p11 <- result2$plots$barplot
p11_path <- file.path(sa_path, "GO barplot miR26a.png")
ggsave(p11_path, p11, width = 8, height = 6, dpi = 600)

p12 <- result2$plots$cnetplot
p12_path <- file.path(sa_path, "GO cnetplot miR26a.png")
ggsave(p12_path, p12, width = 8, height = 6, dpi = 600)

# UTR length v.s. miRNA count
p13 <- ggplot(gene_stats, aes(x = utr_length, y = target.number)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "firebrick") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") + 
  labs(x = "3' UTR length (nt)", 
       y = "Number of miRNA targets",
       title = "Correlation between UTR length and miRNA target count (Liver)") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
p13_path <- file.path(sa_path, "cor UTR miRNA.png")
ggsave(p13_path, p13, width = 8, height = 6, dpi = 600)