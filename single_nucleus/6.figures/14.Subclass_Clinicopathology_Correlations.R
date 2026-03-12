# Clinical/Pathology Correlations -----------------------------------------

library(Seurat)
library(patchwork)
library(dplyr)
library(stringr)
library(Matrix)
library(ggplot2)
library(BPCells)
library("RColorBrewer")
library(gplots)
library(dendextend)

options(future.globals.maxSize = 1e9)

setwd("~/Projects/Human_Kidney/Atlas_V2")
options(Seurat.object.assay.version = "v5")

###Pathology score groupings - Interstitial fibrosis groups
KB <- readRDS("~/datasets/altos-lab-restricted-project-kpmp/blake_LTS/Atlas_V2/scratch/HKAv2_Seurat_01162026.rds")
KB <- subset(KB, tissue_type %in% "Biopsy")
KB <- subset(KB, condition_level1 %in% c("HRT","AKI","CKD"))
KB

#Add pathology Scores
path.meta <- read.delim("KPMP_TIV_Descriptor_Scores_8-29-2024_subset.txt")
path.meta <- path.meta[path.meta$percent_cortex >= 25, ]
colnames(path.meta)
KB <- subset(KB, patient %in% path.meta$patient)

meta <- KB@meta.data
emc <- c("percent_int_fibrosis","percent_tub_atrophy","percent_tub_injury","percent_ab_lum_morph",
         "percent_acellular_casts","percent_int_WBCs","arteriosclerosis","arteriolar_hyalinosis")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- path.meta[,emc[i]][match(meta$patient, path.meta$patient)]
}
colnames(meta)
KB@meta.data <- meta
KB$eGFR <- gsub(" ml/min/1.73m2","",KB$baseline_eGFR_binned)
KB$fibrosis_score <- KB$percent_int_fibrosis

length(unique(KB@meta.data[KB@meta.data$condition_level1 == "HRT",]$patient)) 
#0
length(unique(KB@meta.data[KB@meta.data$condition_level1 == "AKI",]$patient))
#27
length(unique(KB@meta.data[KB@meta.data$condition_level1 == "CKD",]$patient))
#54

## Function to process clinical metadata
process_clinical_metadata <- function(seurat_obj) {
  
  # ===== Process eGFR =====
  convert_egfr <- function(x) {
    x_clean <- str_trim(as.character(x))
    numbers <- as.numeric(str_extract_all(x_clean, "\\d+")[[1]])
    
    if (length(numbers) == 0) return(NA_real_)
    if (length(numbers) == 1) return(numbers[1] + 5)  # For ">60"
    if (length(numbers) == 2) return(mean(numbers))   # For ranges
    return(NA_real_)
  }
  
  assign_ckd_stage <- function(egfr_value) {
    case_when(
      is.na(egfr_value) ~ NA_character_,
      egfr_value >= 90 ~ "CKD1 (Normal)",
      egfr_value >= 60 ~ "CKD2 (Mild)",
      egfr_value >= 45 ~ "CKD3a (Mild-Moderate)",
      egfr_value >= 30 ~ "CKD3b (Moderate-Severe)",
      egfr_value >= 15 ~ "CKD4 (Severe)",
      TRUE ~ "CKD5 (Failure)"
    )
  }
  
  # ===== Process Fibrosis =====
  convert_fibrosis <- function(x) {
    x_clean <- str_trim(as.character(x))
    
    # Handle ranges
    if (str_detect(x_clean, "to|-")) {
      numbers <- as.numeric(str_extract_all(x_clean, "\\d+")[[1]])
      if (length(numbers) == 2) return(mean(numbers))
    }
    
    # Convert single values
    as.numeric(x_clean)
  }
  
  categorize_fibrosis <- function(fib_value) {
    case_when(
      is.na(fib_value) ~ NA_character_,
      fib_value == 0 ~ "None",
      fib_value < 10 ~ "Minimal",
      fib_value < 25 ~ "Mild",
      fib_value < 50 ~ "Moderate",
      TRUE ~ "Severe"
    )
  }
  
  # ===== Apply conversions =====
  message("Converting eGFR...")
  seurat_obj$eGFR_numeric <- sapply(seurat_obj$eGFR, convert_egfr)
  seurat_obj$CKD_stage <- sapply(seurat_obj$eGFR_numeric, assign_ckd_stage)
  
  message("Converting fibrosis scores...")
  seurat_obj$fibrosis_numeric <- sapply(seurat_obj$fibrosis_score, convert_fibrosis)
  seurat_obj$fibrosis_category <- sapply(seurat_obj$fibrosis_numeric, categorize_fibrosis)
  
  # Create ordered factors
  seurat_obj$CKD_stage <- factor(
    seurat_obj$CKD_stage,
    levels = c("CKD1 (Normal)", "CKD2 (Mild)", "CKD3a (Mild-Moderate)", 
               "CKD3b (Moderate-Severe)", "CKD4 (Severe)", "CKD5 (Failure)"),
    ordered = TRUE
  )
  
  seurat_obj$fibrosis_category <- factor(
    seurat_obj$fibrosis_category,
    levels = c("None", "Minimal", "Mild", "Moderate", "Severe"),
    ordered = TRUE
  )
  
  # Summary statistics
  message("\n=== Conversion Summary ===")
  message("\neGFR Distribution:")
  print(summary(seurat_obj$eGFR_numeric))
  message("\nCKD Stage Distribution:")
  print(table(seurat_obj$CKD_stage, useNA = "always"))
  
  message("\nFibrosis Score Distribution:")
  print(summary(seurat_obj$fibrosis_numeric))
  message("\nFibrosis Category Distribution:")
  print(table(seurat_obj$fibrosis_category, useNA = "always"))
  
  # Check for problematic conversions
  n_egfr_na <- sum(is.na(seurat_obj$eGFR_numeric))
  n_fib_na <- sum(is.na(seurat_obj$fibrosis_numeric))
  
  if (n_egfr_na > 0) {
    message(sprintf("\nWarning: %d eGFR values failed to convert", n_egfr_na))
  }
  if (n_fib_na > 0) {
    message(sprintf("Warning: %d fibrosis values failed to convert", n_fib_na))
  }
  
  return(seurat_obj)
}

# Process clinical metadata
KB <- process_clinical_metadata(KB)


###Verify Processing
# Create a verification dataframe
verification <- KB@meta.data %>%
  select(patient, eGFR, eGFR_numeric, CKD_stage, 
         fibrosis_score, fibrosis_numeric, fibrosis_category) %>%
  distinct()

# View unique combinations
verification %>%
  arrange(eGFR_numeric) %>%
  distinct(eGFR, eGFR_numeric, CKD_stage)

verification %>%
  arrange(fibrosis_numeric) %>%
  distinct(fibrosis_score, fibrosis_numeric, fibrosis_category)

# Check for any failed conversions
verification %>%
  filter(is.na(eGFR_numeric) | is.na(fibrosis_numeric))


###Correlation Analyses
library(broom)

# Calculate cell type proportions
cell_props <- KB@meta.data %>%
  group_by(patient, v2.subclass.l3) %>%
  summarise(n_cells = n(), .groups = 'drop') %>%
  group_by(patient) %>%
  mutate(proportion = n_cells / sum(n_cells))

# Add numeric clinical data
clinical_data <- KB@meta.data %>%
  select(patient, eGFR_numeric, fibrosis_numeric, 
         CKD_stage, fibrosis_category) %>%
  distinct()

cell_props <- left_join(cell_props, clinical_data, by = "patient")


# Check sample sizes per cell type
sample_sizes <- cell_props %>%
  group_by(v2.subclass.l3) %>%
  summarise(
    n_patients = n(),
    n_complete_eGFR = sum(!is.na(eGFR_numeric) & !is.na(proportion)),
    n_complete_fibrosis = sum(!is.na(fibrosis_numeric) & !is.na(proportion))
  )

print(sample_sizes, n = 128)


###Run analyses for numeric values
results <- cell_props %>%
  group_by(v2.subclass.l3) %>%
  summarise(
    # Correlations (continuous)
    cor_eGFR = cor(proportion, eGFR_numeric, method = "spearman", use = "complete.obs"),
    cor_fibrosis = cor(proportion, fibrosis_numeric, method = "spearman", use = "complete.obs"),
    
    # Linear regression (continuous)
    slope_eGFR = coef(lm(proportion ~ eGFR_numeric))[2],
    pval_eGFR = summary(lm(proportion ~ eGFR_numeric))$coefficients[2, 4],
    r2_eGFR = summary(lm(proportion ~ eGFR_numeric))$r.squared,
    
    slope_fibrosis = coef(lm(proportion ~ fibrosis_numeric))[2],
    pval_fibrosis = summary(lm(proportion ~ fibrosis_numeric))$coefficients[2, 4],
    r2_fibrosis = summary(lm(proportion ~ fibrosis_numeric))$r.squared,
    
    # Joint model
    model_both = list(broom::tidy(lm(proportion ~ eGFR_numeric + fibrosis_numeric))),
    
    n_patients = n()
  ) %>%
  mutate(
    padj_eGFR = p.adjust(pval_eGFR, method = "BH"),
    padj_fibrosis = p.adjust(pval_fibrosis, method = "BH")
  )

saveRDS(results, file = "Clinical_Correlation_Results_CellTypeProportions_SubclassL3_Biopsy_eGFR_Fibrosis_01072026.rds")
results.tab <- as.data.frame(results)
results.tab <- results.tab %>%
  select(-model_both)
write.csv(results.tab, file = "Clinical_Correlation_Results_CellTypeProportions_SubclassL3_Biopsy_eGFR_Fibrosis_01072026.csv")

# View results
results %>%
  select(v2.subclass.l3, cor_eGFR, padj_eGFR, cor_fibrosis, padj_fibrosis, n_patients) %>%
  arrange(padj_eGFR)%>%
  print(n = Inf)

summary(results$cor_eGFR)
summary(results$cor_fibrosis)

results %>%
  select(v2.subclass.l3, cor_eGFR, pval_eGFR, pval_eGFR, cor_fibrosis, pval_fibrosis) %>%
  arrange(pval_fibrosis) %>%
  print(n = Inf)


##Heatmap for all cell types
# Prepare data for heatmap
cor_matrix <- results %>%
  select(cor_eGFR, cor_fibrosis) 

cor_matrix <- as.data.frame(cor_matrix)
colnames(cor_matrix) <- c("eGFR", "Fibrosis Score")
rownames(cor_matrix) <- results$v2.subclass.l3

# Create significance matrix for annotations
sig_matrix <- results %>%
  select(pval_eGFR, pval_fibrosis)
sig_matrix <- as.data.frame(sig_matrix)
colnames(sig_matrix) <- c("eGFR", "Fibrosis Score")
rownames(sig_matrix) <- results$v2.subclass.l3

order <- c("POD","dPOD","PEC","PT-S1","dPT-S1","PT-S2","dPT-S2","PT-S3","dPT-S3",
           "aPT2","aPT1","aPT-S1/S2","frPT-S1/S2","aPT-S3","frPT-S3","cycPT","dPT","DTL2","aDTL2","DTL1","DTL3",
           "dDTL","ATL","dATL","M-TAL","dM-TAL","C/M-TAL-A","C-TAL-A","MD","C/M-TAL-B","C-TAL-B","dC-TAL",
           "aTAL1","aTAL2","frTAL","cycTAL","DCT1","DCT2","dDCT","frDCT","aDCT","CNT","dCNT","aCNT",
           "CNT-PC","dCNT-PC","CCD-PC","OMCD-PC","dOMCD-PC","IMCD","dIMCD","PapE","CCD-IC-A","dCCD-IC-A",
           "OMCD-IC-A","dOMCD-IC-A","tOMCD-PC-IC","IC-B","EC-GC","dEC-GC","aEC-GC","EC-GC-FILIP1+","EC-AA",
           "EC-DVR","M-EC-PTC","EC-AVR","dEC-AVR","iaEC-AVR","EC-V","EC-PCV","C-EC-PTC","angEC-PTC",
           "EC-EA","dEC-PTC","infEC-PTC","iaEC-PTC","EC-LYM","cycEC","IM-FIB","OM-FIB","dIM-FIB","dOM-FIB",
           "C/M-FIB","IM-pvMYOF","C-FIB","C-FIB-PATH","C-FIB-OSMRlo","C-FIB-OSMRhi","C-MYOF","dFIB",
           "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF","MC","REN","VSMC","M-VSMC/P","VSMC/P","dVSMC",
           "Ad","B","PL","Naive Th","MAIT","ILC3","T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY",
           "MAST","resMAC-LYVE1+","resMAC-HLAIIhi","MON","moMAC-HBEGF+","moMAC-CXCL10+","moFAM","moMAC-C3+",
           "cDC2","ncMON","mDC","cDC1","pDC","N","cycT","cycMAC","SC/NEU")


pdf(file='QC_Plots/Composition_Correlation_Analysis_Heatmap_ordered_subclass_level3_biopsy.pdf',width=30,height=6)
heatmap.2(as.matrix(t(cor_matrix[order,])),col=brewer.pal(11,"RdBu"),scale="none", trace="none", 
          density.info= "none", cexRow = 1, margins = c(8, 8), Rowv = NA, Colv = NA,
          #ColSideColors = ColSideColors,
          cellnote = ifelse(as.matrix(t(sig_matrix[order,])) < 0.05, "*", ""),
          notecol = 'black', notecex = 2)
dev.off() 



pdf(file='QC_Plots/Composition_Analysis_t-statistic_Heatmap_Clinical_ordered_subclass_level3_biopsy_B.pdf',width=30,height=3)
heatmap.2(as.matrix(t(cor_matrix[order,])),col=brewer.pal(11,"RdBu"),scale="none", trace="none", 
          density.info= "none", cexRow = 1, margins = c(8, 8), Rowv = NA, Colv = NA,
          #ColSideColors = ColSideColors,
          cellnote = ifelse(as.matrix(t(sig_matrix[order,])) < 0.05, "*", ""),
          notecol = 'black', notecex = 2)
dev.off() 



##Plot select cell types

# Function to plot with proper labels
plot_clinical_correlation <- function(data, ct, variable, var_label) {
  subset_data <- data %>% filter(v2.subclass.l3 == ct)
  
  # Calculate stats
  cor_val <- cor(subset_data$proportion, subset_data[[variable]], 
                 method = "spearman", use = "complete.obs")
  lm_fit <- lm(as.formula(paste("proportion ~", variable)), data = subset_data)
  
  ggplot(subset_data, aes(x = .data[[variable]], y = proportion)) +
    geom_point(size = 3, alpha = 0.6, color = "steelblue") +
    geom_smooth(method = "lm", se = TRUE, color = "red", fill = "pink") +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("r = %.3f\nSlope = %.4f\nR² = %.3f\np = %.3g",
                             cor_val,
                             coef(lm_fit)[2],
                             summary(lm_fit)$r.squared,
                             summary(lm_fit)$coefficients[2, 4]),
             hjust = 1.1, vjust = 1.1, size = 3.5) +
    labs(title = paste(ct, "vs", var_label),
         x = var_label,
         y = "Cell Type Proportion") +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", size = 12))
}

# frTAL vs eGFR and Fibrosis Score
p1 <- plot_clinical_correlation(cell_props, "frTAL", 
                                "eGFR_numeric", "eGFR (ml/min/1.73m²)")
p2 <- plot_clinical_correlation(cell_props, "frTAL", 
                                "fibrosis_numeric", "Fibrosis Score (%)")

pdf(file='QC_Plots/Composition_Analysis_frTAL_vs_eGFR_and_FibrosisScore_biopsy.pdf',width=8,height=4)
p1 + p2
dev.off()



# frPT-S1/S2 vs eGFR and Fibrosis Score
p1 <- plot_clinical_correlation(cell_props, "frPT-S1/S2", 
                                "eGFR_numeric", "eGFR (ml/min/1.73m²)")
p2 <- plot_clinical_correlation(cell_props, "frPT-S1/S2", 
                                "fibrosis_numeric", "Fibrosis Score (%)")

pdf(file='QC_Plots/Composition_Analysis_frPT-S1S2_vs_eGFR_and_FibrosisScore_biopsy.pdf',width=8,height=4)
p1 + p2
dev.off()



# frPT-S3 vs eGFR and Fibrosis Score
p1 <- plot_clinical_correlation(cell_props, "frPT-S3", 
                                "eGFR_numeric", "eGFR (ml/min/1.73m²)")
p2 <- plot_clinical_correlation(cell_props, "frPT-S3", 
                                "fibrosis_numeric", "Fibrosis Score (%)")

pdf(file='QC_Plots/Composition_Analysis_frPT-S3_vs_eGFR_and_FibrosisScore_biopsy.pdf',width=8,height=4)
p1 + p2
dev.off()


