# Clear workspace
library(dplyr)
library(SomaDataIO)
library(purrr)
library(tidyverse)
library(readxl)
library(writexl)
library(httr)

# The first script. Data cleaning
setwd("~/Desktop/有事/Harvard/Research/TracyYoung-Pearse_Proteomics/Zach/code")
source("./functions.R")

################### Data Cleaning ###################
# Load data
iA_origin <- read_excel("../data/BAG3 iA Raw Data.xlsx", sheet = "updated")
# iN_df <- read_excel("../Zach/data/BAG3 iN Raw Proteomic Data.xlsx")

# Define the mapping for column renaming
column_mapping <- c(
  "127N" = "BR24_WT_r1", "127C" = "BR24_WT_r2", "128N" = "BR24_WT_r3", "128C" = "BR24_WT_r4",
  "129N" = "BR24_KO_r1", "129C" = "BR24_KO_r2", "130N" = "BR24_KO_r3", "130C" = "BR24_KO_r4",
  "131N" = "BR33_WT_r1", "131C" = "BR33_WT_r2", "132N" = "BR33_WT_r3", "132C" = "BR33_WT_r4",
  "133N" = "BR33_KO_r1", "133C" = "BR33_KO_r2", "134N" = "BR33_KO_r3", "134C" = "BR33_KO_r4"
)

# Select and rename columns
iA_df <- iA_origin %>%
  dplyr::select(Accession, Description, starts_with("Abundances (Normalized): F3:")) %>%
  rename_with(~ column_mapping[sub(".*: ([^,]+).*", "\\1", .)], 
              starts_with("Abundances (Normalized): F3:")) %>%
  mutate(Target = str_extract(Description, "^.*?(?=\\sOS=)"),
         EntrezGeneSymbol = str_extract(Description, "(?<=GN=)\\w+")) %>%
  dplyr::select(Accession, Target, EntrezGeneSymbol, everything())

prot_cand <- iA_df %>% dplyr::select(Accession, Target, EntrezGeneSymbol)
saveRDS(prot_cand, file = "../data/prot_cand.rds")

# Transpose data frame
sample_info <- tibble(
  sample = c("BR24_WT_r1", "BR24_WT_r2", "BR24_WT_r3", "BR24_WT_r4",
             "BR24_KO_r1", "BR24_KO_r2", "BR24_KO_r3", "BR24_KO_r4",
             "BR33_WT_r1", "BR33_WT_r2", "BR33_WT_r3", "BR33_WT_r4",
             "BR33_KO_r1", "BR33_KO_r2", "BR33_KO_r3", "BR33_KO_r4"),
  BRID = rep(c("BR24", "BR33"), each = 8),
  Condition = rep(rep(c("WT", "KO"), each = 4), 2),
  Replicate = rep(1:4, 4))

transposed_df <- iA_df %>%
  dplyr::select(-Target, -EntrezGeneSymbol, -Description) %>% 
  pivot_longer(cols = -Accession, names_to = "sample", values_to = "value") %>%
  pivot_wider(names_from = Accession, values_from = value)

iA_df_norm <- transposed_df %>%
  left_join(sample_info, by = "sample") %>%
  dplyr::select(sample, BRID, Condition, Replicate, everything()) %>%
  arrange(BRID, Condition, Replicate)

saveRDS(iA_df_norm, file = "../data/norm_iA.rds")
iA_df_norm <- readRDS("../data/norm_iA.rds")

################### Data Preparation ###################

pdf("../result/normalization.pdf", width = 15, height = 6)

# Compare normalization
before_normalization <- iA_origin %>% 
  dplyr::select(Accession, starts_with("Abundance: F3:")) %>%
  rename_with(~ column_mapping[sub(".*: ([^,]+).*", "\\1", .)], starts_with("Abundance: F3:")) %>%
  pivot_longer(cols = -Accession, names_to = "sample", values_to = "value") %>%
  pivot_wider(names_from = Accession, values_from = value) %>%
  dplyr::select(-sample)
protein_mean <- apply(before_normalization, 2, mean, na.rm = TRUE)
protein_sd <- apply(before_normalization, 2, sd, na.rm = TRUE)

after_normalization <- iA_df_norm %>% dplyr::select(-sample, -BRID, -Condition, -Replicate)
clean_mean <- apply(after_normalization, 2, mean, na.rm = TRUE)
clean_sd <- apply(after_normalization, 2, sd, na.rm = TRUE)

cleanData <- iA_df_norm %>%
  mutate(across(-c(sample, BRID, Condition, Replicate), ~log10(. + 1))) %>%  # log10-transform all except specified columns
  mutate(across(-c(sample, BRID, Condition, Replicate), cs))         # center/scale all except specified columns

cs_normalization <- cleanData %>% dplyr::select(-sample, -BRID, -Condition, -Replicate)
cs_mean <- apply(cs_normalization, 2, mean, na.rm = TRUE)
cs_sd <- apply(cs_normalization, 2, sd, na.rm = TRUE)

par(mfrow = c(1, 3))  # Create a 1x2 plot layout
plot(protein_mean, protein_sd, main = "Before Normalization", xlab = "Mean", ylab = "Standard Deviation", pch = 16)
plot(clean_mean, clean_sd, main = "After Normalization", xlab = "Mean", ylab = "Standard Deviation", pch = 16)
plot(cs_mean, cs_sd, main = "After Center/Scale", xlab = "Mean", ylab = "Standard Deviation", pch = 16)

na_columns <- names(cleanData)[sapply(cleanData, function(x) all(is.na(x)))]
cat("Number of Columns with all NA values:\n")
print(length(na_columns))

cleanData <- cleanData[, !sapply(cleanData, function(x) all(is.na(x)))]
# saveRDS(cleanData, "../data/cleanData.rds")
dev.off()

pdf("../result/pca.pdf", width = 6, height = 6)
# 2. PCA plot for BR columns
input <- cleanData
strata1 <- "Condition" # Treatment or Control
strata2 <- "BRID" # BR24 or BR33
strata1_num <- length(unique(input[[strata1]]))
strata2_num <- length(unique(input[[strata2]]))

metadata <- input %>% dplyr::select("sample", "Condition", "BRID", "Replicate")
protein_columns <- setdiff(colnames(input), c("sample", "Condition", "BRID", "Replicate"))
input <- input %>% 
  dplyr::select(sample, all_of(protein_columns)) %>%
  column_to_rownames("sample")

pca_result <- prcomp(input, center = FALSE, scale. = FALSE)

# Calculate the percentage of variance explained by PC1 and PC2
variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
pc1_variance <- round(variance_explained[1] * 100, 2)
pc2_variance <- round(variance_explained[2] * 100, 2)
pc3_variance <- round(variance_explained[3] * 100, 2)

pca_data <- as.data.frame(pca_result$x)
pca_data$sample <- rownames(pca_data)
pca_data <- pca_data %>% left_join(metadata, by = "sample")
pca_data[[strata1]] <- factor(pca_data[[strata1]])
pca_data[[strata2]] <- factor(pca_data[[strata2]])

# Define colors and shapes
strata2_colors <- rainbow_hcl(strata2_num)
strata1_shapes <- c(16, 17, 18, 15, 3, 8, 9, 10, 11, 12)[1:strata1_num]

ggplot(pca_data, aes(x = PC1, y = PC2, color = .data[[strata2]], shape = .data[[strata1]])) +
  geom_point(size = 8, alpha = 0.6) +
  ggtitle("PCA of Protein Data (PC1 vs PC2)") +
  xlab(paste0("Principal Component 1 (", pc1_variance, "%)")) +
  ylab(paste0("Principal Component 2 (", pc2_variance, "%)")) +
  scale_color_manual(values = strata2_colors, name = strata2) +
  scale_shape_manual(values = strata1_shapes, name = strata1) +
  theme_minimal()

# Different colors for strata2 (timepoints) and different shapes for strata1 (group) values
ggplot(pca_data, aes(x = PC3, y = PC2, color = .data[[strata2]], shape = .data[[strata1]])) +
  geom_point(size = 8, alpha = 0.6) +
  ggtitle("PCA of Protein Data (PC3 vs PC2)") +
  xlab(paste0("Principal Component 3 (", pc3_variance, "%)")) +
  ylab(paste0("Principal Component 2 (", pc2_variance, "%)")) +
  scale_color_manual(values = strata2_colors, name = strata2) +
  scale_shape_manual(values = strata1_shapes, name = strata1) +
  theme_minimal()
dev.off()

# If see Sample Normalization (non-transposed) 
# before_normalization <- iA_df %>% 
#   dplyr::select(Accession, Description, starts_with("Abundance: F3:")) %>%
#   rename_with(~ column_mapping[sub(".*: ([^,]+).*", "\\1", .)], 
#               starts_with("Abundance: F3:")) %>%
#   dplyr::select(-Accession, -Description)
# 
# after_normalization <- iA_df_norm %>%
#   dplyr::select(-Accession, -Target, -EntrezGeneSymbol)
# 
# iA_df_cs <- iA_df_norm %>%
#   mutate(across(-c(Accession, Target, EntrezGeneSymbol), ~log10(. + 1))) %>%  # log10-transform all except specified columns
#   mutate(across(-c(Accession, Target, EntrezGeneSymbol), cs))         # center/scale all except specified columns
# 
# after_cs <- iA_df_cs %>%
#   dplyr::select(-Accession, -Target, -EntrezGeneSymbol)
# 
# par(mfrow = c(1, 3))  # Create a 1x3 plot layout
# 
# protein_mean <- apply(before_normalization, 2, mean, na.rm = TRUE)
# protein_sd <- apply(before_normalization, 2, sd, na.rm = TRUE)
# plot(protein_mean, protein_sd, main = "Before Normalization", xlab = "Mean", ylab = "Standard Deviation", pch = 16)
# 
# clean_mean <- apply(after_normalization, 2, mean, na.rm = TRUE)
# clean_sd <- apply(after_normalization, 2, sd, na.rm = TRUE)
# plot(clean_mean, clean_sd, main = "After Normalization", xlab = "Mean", ylab = "Standard Deviation", pch = 16)
# 
# cs_mean <- apply(after_cs, 2, mean, na.rm = TRUE)
# cs_sd <- apply(after_cs, 2, sd, na.rm = TRUE)
# plot(cs_mean, cs_sd, main = "After Center/Scale", xlab = "Mean", ylab = "Standard Deviation", 
#      pch = 16, xlim = c(-5e-15, 5e-15), ylim = c(0.9, 1.1))

################################ Analysis ######################################
# Split data to two groups
cleanData <- readRDS("../data/cleanData.rds")
prot_cand <- readRDS("../data/prot_cand.rds")
source("./functions.R")

# Define analysis function
perform_analysis <- function(data, brid) {
  reshaped_data <- data %>%
    as.data.frame() %>% 
    dplyr::select(-BRID, -Condition, -Replicate) %>%
    pivot_longer(cols = -sample, names_to = "protein", values_to = "expression") %>%
    pivot_wider(names_from = sample, values_from = expression) %>%
    dplyr::select(protein, everything())
  
  method.list <- c("Limma", "ROTS")
  
  plotp <- paste0("../result/plots_", brid, ".pdf")
  pdf(plotp, width = 6, height = 6)
  
  # de_results <- lapply(method.list, function(m) {
  #   source(paste("Methods/", m, ".R", sep = ""))
  #   print(m)
  #   set.seed(1)
  #   getDE(reshaped_data)
  # })
  # names(de_results) <- method.list
  # saveRDS(de_results, file = paste0("../data/de_results_", brid, ".rds"))
  de_results <- readRDS(paste0("../data/de_results_", brid, ".rds"))
    
  prot_list <- lapply(method.list, function(m) {
    de_res <- de_results[[m]]
    prot <- add_EntrezGeneSymbol(de_res, prot_cand)
    write_xlsx(prot, paste0("../data/", m, "_res_", brid, ".xlsx"))
    prot
  })
  names(prot_list) <- method.list
  
  upset_result <- create_upset_plot(method.list, brid)
  overlap_df <- upset_result$overlapping_proteins
  
  grid.draw(upset_result$venn_plot)
  print(upset_result$upset_plot)
  
  write_xlsx(overlap_df, paste0("../result/limmaROTS_overlap_proteins_", brid, ".xlsx"))
  
  print(perform_go_analysis_and_plot(overlap_df, plotp, method.list))
  dev.off()
  
  return(overlap_df)
}

iA_br24 <- cleanData %>% filter(str_detect(sample, "^BR24"))
iA_br33 <- cleanData %>% filter(str_detect(sample, "^BR33"))
results_br24 <- perform_analysis(iA_br24, "BR24")
results_br33 <- perform_analysis(iA_br33, "BR33")

br24_33 <- list(br24=results_br24, br33=results_br33)

# Prepare Upset plot dataset
br2433_clean <- lapply(br24_33, function(x) x$EntrezGeneSymbol[!is.na(x$EntrezGeneSymbol)])
mlist <- UpSetR::fromList(br2433_clean)

# Create UpSet plot
upset_plot <- upset(mlist, nsets = 2, 
                    text.scale = c(2, 2, 1.5, 1.5, 2, 2),
                    line.size = 1.5,
                    nintersects = NA, 
                    order.by = "freq")

# Create Venn diagram
# color_palette <- c("#F7FEF0", "#E9F7E5", "#CEEFCC", "#BFE8C1", "#BCF4C5", "#92C2A6",
#                    "#D6F6FF", "#ACEEFE", "#6FC8CA", "#58B8D1", "#3492B2", "#04579B")
color_palette <- c("#F7FEF0", "#BFE8C1", "#58B8D1", "#04579B")
venn_plot <- venn.diagram(
  x = br2433_clean,
  filename = NULL,
  fill = color_palette[1:2],
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  fontfamily = "sans",
  disable.logging = TRUE
)

overlapping_proteins <- as.data.frame(intersect(results_br24$EntrezGeneSymbol, results_br33$EntrezGeneSymbol))
write_xlsx(overlapping_proteins, "../result/overlap_br24_33.xlsx")

significant_genes <- bitr(overlapping_proteins$`intersect(results_br24$EntrezGeneSymbol, results_br33$EntrezGeneSymbol)`, 
                          fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go_result <- enrichGO(gene = significant_genes$ENTREZID,
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.5,
                      qvalueCutoff = 0.5)
plot <- dotplot(go_result, showCategory = 15) + 
  ggtitle(paste0("Top 15 GO BP for BR24 and BR33 Overlap"))

pdf("../result/br24_33_plots.pdf", width = 8, height = 8)
grid.draw(venn_plot)
print(upset_plot)
print(plot)
dev.off()
############################### More Plots #####################################
# Volcano Plot
prot_cand <- readRDS("../data/prot_cand.rds")
source("./functions.R")

method <- "ROTS"
type <- "iA"

volcano_plot <- function(brid) {
  prot_list <- read_xlsx(paste0("../data/", method, "_res_", brid, ".xlsx"))

  # Use all proteins as prot_interest
  prot_interest <- prot_list %>%
    group_by(EntrezGeneSymbol) %>%
    filter(fdr == min(fdr)) %>%
    ungroup() %>%
    pull(ID)

  # Get fold change matrices for all combinations
  fold_change_matrices <- get_fold_change_matrices(type, brid, prot_cand, prot_interest)

  # Create a data frame for the volcano plot
  volcano_data <- data.frame(
    EntrezGeneSymbol = fold_change_matrices$EntrezGeneSymbol,
    log2FoldChange = fold_change_matrices$fold_change,
    padj = prot_list$fdr[match(fold_change_matrices$EntrezGeneSymbol, prot_list$EntrezGeneSymbol)]
  )
  volcano_plot <- volcano_fun(type, brid, volcano_data)
  
  return(volcano_plot)
}

# Perform analysis for both BRIDs
plot_BR24 <- volcano_plot("BR24")
plot_BR33 <- volcano_plot("BR33")

plotp <- paste0("../result/", type, "_volcano.pdf")
pdf(plotp, width = 6, height = 8)
print(plot_BR24)
print(plot_BR33)
dev.off()



