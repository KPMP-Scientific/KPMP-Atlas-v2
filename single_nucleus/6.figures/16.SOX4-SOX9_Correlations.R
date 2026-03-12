# PT-S1-2 FR Trajectory Sox4/Sox9 Activity Correlations ---------------------------------
library(Seurat)
library(Signac)
library(patchwork)
library(dplyr)
library(Matrix)
library(ggplot2)
library(BPCells)

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

setwd("~/Projects/Human_Kidney/Atlas_V2/")

load("~/datasets/altos-lab-restricted-project-kpmp/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_aPT-S1-2_Trajectory_GRN_0424-newData.RDA")
tf_target_connections <- df.grn[, c("tf", "gene", "correlation", "p_value")]
colnames(tf_target_connections) <- c("TF", "target_gene", "weight", "p_value")

#/home/blake/hsKidAt/blake_LTS/Atlas_V2/scratch/full_kidney_count_set_0424
counts <- open_matrix_dir(dir = "/home/blake/datasets/altos-lab-restricted-project-kpmp/blake_LTS/Atlas_V2/scratch/full_kidney_count_set_0424")
obj[["RNA"]]$counts <- counts
obj <- NormalizeData(obj, assay = "RNA")
obj[["RNA"]]$counts <- as(obj[["RNA"]]$counts, Class = "dgCMatrix")
seurat_obj <- obj

# Extract SOX4 and SOX9 targets
sox4_targets <- tf_target_connections %>%
  filter(TF == "SOX4", abs(weight) > 0.8) %>%
  pull(target_gene) %>%
  unique()

sox9_targets <- tf_target_connections %>%
  filter(TF == "SOX9", abs(weight) > 0.8) %>%
  pull(target_gene) %>%
  unique()

cat(sprintf("SOX4 targets: %d\n", length(sox4_targets)))
cat(sprintf("SOX9 targets: %d\n", length(sox9_targets)))
cat(sprintf("Shared targets: %d\n", length(intersect(sox4_targets, sox9_targets))))
obj

###Extract Expression, Activity, and Pseudotime Data
# Function to extract all relevant data
extract_tf_data <- function(seurat_obj, tf_name, targets, 
                            pseudotime_col = "Trajectory",
                            activity_assay = "chromvar") {
  
  # Get pseudotime
  pseudotime <- seurat_obj@meta.data[[pseudotime_col]]
  
  # Get TF expression (RNA)
  if (tf_name %in% rownames(seurat_obj@assays$RNA)) {
    tf_expression <- GetAssayData(seurat_obj, assay = "RNA", 
                                  layer = "data")[tf_name, ]
    tf_expression <- as(tf_expression, "dgCMatrix")
  } else {
    warning(sprintf("%s not found in RNA assay", tf_name))
    tf_expression <- rep(NA, ncol(seurat_obj))
  }
  
  # Get TF activity (chromVAR or motif activity)
  # Adjust based on your scMEGA output structure
  if (activity_assay %in% names(seurat_obj@assays)) {
    # Try to find the motif/TF in the activity matrix
    motif_names <- rownames(seurat_obj@assays[[activity_assay]])
    names(motif_names) <- ConvertMotifID(object = seurat_obj, id = motif_names)
    matching_motif <- as.character(motif_names[tf_name])
    
    if (length(matching_motif) > 0) {
      tf_activity <- GetAssayData(seurat_obj, assay = activity_assay, 
                                  layer = "data")[matching_motif[1], ]
    } else {
      warning(sprintf("No activity score found for %s", tf_name))
      tf_activity <- rep(NA, ncol(seurat_obj))
    }
  } else {
    warning(sprintf("Activity assay '%s' not found", activity_assay))
    tf_activity <- rep(NA, ncol(seurat_obj))
  }
  
  # Get target gene expression (average across targets)
  targets_in_data <- intersect(targets, rownames(seurat_obj@assays$RNA))
  
  if (length(targets_in_data) > 0) {
    target_expr <- GetAssayData(seurat_obj, assay = "RNA", 
                                layer = "data")[targets_in_data, , drop = FALSE]
    target_expr <- as(target_expr, "dgCMatrix")
    avg_target_expression <- colMeans(as.matrix(target_expr))
  } else {
    warning(sprintf("No targets found for %s", tf_name))
    avg_target_expression <- rep(NA, ncol(seurat_obj))
  }
  
  # Combine into dataframe
  df <- data.frame(
    cell = colnames(seurat_obj),
    pseudotime = pseudotime,
    tf_expression = as.numeric(tf_expression),
    tf_activity = as.numeric(tf_activity),
    avg_target_expression = avg_target_expression,
    stringsAsFactors = FALSE
  )
  
  return(df)
}

# Extract data for both TFs
sox4_data <- extract_tf_data(obj, "SOX4", sox4_targets)
sox9_data <- extract_tf_data(obj, "SOX9", sox9_targets)

# Merge the dataframes
combined_data <- sox4_data %>%
  select(cell, pseudotime, 
         SOX4_expression = tf_expression,
         SOX4_activity = tf_activity,
         SOX4_target_expr = avg_target_expression) %>%
  left_join(
    sox9_data %>%
      select(cell,
             SOX9_expression = tf_expression,
             SOX9_activity = tf_activity,
             SOX9_target_expr = avg_target_expression),
    by = "cell"
  )

head(combined_data)

###Correlation Analysis
# Calculate all pairwise correlations
calculate_correlations <- function(data, method = "spearman") {
  
  # Select numeric columns for correlation
  cor_cols <- c("pseudotime", "SOX4_expression", "SOX4_activity", 
                "SOX4_target_expr", "SOX9_expression", "SOX9_activity", 
                "SOX9_target_expr")
  
  cor_data <- data %>%
    select(all_of(cor_cols)) %>%
    na.omit()
  
  # Calculate correlation matrix
  cor_matrix <- cor(cor_data, method = method)
  
  # Calculate p-values
  n <- nrow(cor_data)
  p_matrix <- matrix(NA, nrow = ncol(cor_data), ncol = ncol(cor_data))
  rownames(p_matrix) <- colnames(p_matrix) <- colnames(cor_data)
  
  for (i in 1:ncol(cor_data)) {
    for (j in 1:ncol(cor_data)) {
      if (i != j) {
        test <- cor.test(cor_data[, i], cor_data[, j], method = method)
        p_matrix[i, j] <- test$p.value
      }
    }
  }
  
  return(list(
    correlation = cor_matrix,
    p_value = p_matrix,
    n_cells = n
  ))
}

# Calculate correlations
correlations <- calculate_correlations(combined_data, method = "spearman")

# Print key correlations
cat("\n=== Key Correlations (Spearman) ===\n")
cat(sprintf("SOX4 vs SOX9 expression: r=%.3f, p=%.3e\n",
            correlations$correlation["SOX4_expression", "SOX9_expression"],
            correlations$p_value["SOX4_expression", "SOX9_expression"]))
cat(sprintf("SOX4 vs SOX9 activity: r=%.3f, p=%.3e\n",
            correlations$correlation["SOX4_activity", "SOX9_activity"],
            correlations$p_value["SOX4_activity", "SOX9_activity"]))
cat(sprintf("SOX4 vs SOX9 target expression: r=%.3f, p=%.3e\n",
            correlations$correlation["SOX4_target_expr", "SOX9_target_expr"],
            correlations$p_value["SOX4_target_expr", "SOX9_target_expr"]))



### Visualize correlation matrix
# Create a nice correlation heatmap with significance stars
library(corrplot)

# Prepare labels with significance stars
add_significance_stars <- function(cor_mat, p_mat) {
  sig_mat <- matrix("", nrow = nrow(p_mat), ncol = ncol(p_mat))
  sig_mat[p_mat < 0.001] <- "***"
  sig_mat[p_mat >= 0.001 & p_mat < 0.01] <- "**"
  sig_mat[p_mat >= 0.01 & p_mat < 0.05] <- "*"
  return(sig_mat)
}

sig_stars <- add_significance_stars(correlations$correlation, 
                                    correlations$p_value)

# Plot correlation heatmap
corrplot(correlations$correlation,
         method = "color",
         type = "upper",
         order = "original",
         addCoef.col = "black",
         tl.col = "black",
         tl.srt = 45,
         diag = FALSE,
         col = colorRampPalette(c("red", "white", "blue"))(200),
         title = "SOX4 and SOX9 Multi-level Correlations",
         mar = c(0, 0, 2, 0))



###Scatter Plots: Expression vs Activity vs Targets
# Create comprehensive scatter plot matrix
library(GGally)

plot_data <- combined_data %>%
  select(SOX4_expression, SOX4_activity, SOX4_target_expr,
         SOX9_expression, SOX9_activity, SOX9_target_expr,
         pseudotime) %>%
  na.omit()

# Rename for cleaner labels
colnames(plot_data) <- c("SOX4\nExpression", "SOX4\nActivity", 
                         "SOX4\nTargets", "SOX9\nExpression", 
                         "SOX9\nActivity", "SOX9\nTargets", 
                         "Pseudotime")

# Individual detailed scatter plots
create_scatter_with_trajectory <- function(data, x_var, y_var, 
                                           x_label, y_label, title) {
  
  cor_test <- cor.test(data[[x_var]], data[[y_var]], 
                       method = "spearman", exact = FALSE)
  
  p <- ggplot(data, aes_string(x = x_var, y = y_var)) +
    geom_point(aes(color = pseudotime), alpha = 0.6, size = 1.5) +
    geom_smooth(method = "lm", color = "red", se = TRUE, alpha = 0.2) +
    scale_color_viridis_c(name = "Pseudotime") +
    labs(
      x = x_label,
      y = y_label,
      title = title,
      subtitle = sprintf("Spearman r = %.3f, p = %.3e", 
                         cor_test$estimate, cor_test$p.value)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
  
  return(p)
}

# Key comparisons
p1 <- create_scatter_with_trajectory(
  combined_data, "SOX4_activity", "SOX9_activity",
  "SOX4 Activity", "SOX9 Activity",
  "SOX4 vs SOX9 Activity"
)

p2 <- create_scatter_with_trajectory(
  combined_data, "SOX4_target_expr", "SOX9_target_expr",
  "SOX4 Target Expression", "SOX9 Target Expression",
  "SOX4 vs SOX9 Target Expression"
)


# Combine plots
pdf(file='regulatory/scMEGA_PT-S1-2_fr-Trajectory_SOX4-SOX9_Activity-TargetExpr_Correlations_Scatterplots.pdf',width=10,height=5)
(p1 + p2)
dev.off()
