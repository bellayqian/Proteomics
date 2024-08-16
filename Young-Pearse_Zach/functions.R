# Load Library
library(ggplot2)
library(dplyr)
library(SomaDataIO)
library(purrr)
library(tidyverse)
library(pheatmap)
library(tidyverse)
library(writexl)
library(readxl)
library(colorspace)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pROC)
library(VennDiagram)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(RolDE)
library(UpSetR)

# Function to add EntrezGeneSymbol to results
add_EntrezGeneSymbol <- function(data, prot_cand){
  library(writexl)
  if ("adj.P.Val" %in% colnames(data)) {
    colnames(data)[colnames(data) == "adj.P.Val"] <- "fdr"
  }
  if ("Feature.ID" %in% colnames(data)) {
    colnames(data)[colnames(data) %in% "Feature.ID"] <- "ID"
  }
  prot_list <- data %>%
    left_join(prot_cand %>% dplyr::select("Accession", Target, EntrezGeneSymbol), 
              by = c("ID" = "Accession")) %>%
    dplyr::select(ID, EntrezGeneSymbol, Target, fdr)
  return(prot_list)
}

# Function to create UpSet plot
create_upset_plot <- function(method.list, brid) {
  m1 <- method.list[1]
  m2 <- method.list[2]
  
  top100_list <- lapply(method.list, function(m) {
    top100 <- read_excel(paste0("../data/", m, "_res_", brid, ".xlsx")) %>% 
      filter(fdr < 0.001) %>%
      filter(!grepl("No protein", Target)) %>%
      filter(!grepl("MOUSE", Target)) %>%
      arrange(fdr)
    top100 <- top100 %>% dplyr::select(ID, EntrezGeneSymbol, Target, fdr)
  })
  names(top100_list) <- method.list
  
  # Create a dataframe with overlapping proteins and their FDR values
  overlapping_proteins <- intersect(top100_list[[m1]]$ID, top100_list[[m2]]$ID)
  result_df <- data.frame(
    ID = overlapping_proteins,
    EntrezGeneSymbol = top100_list[[m1]]$EntrezGeneSymbol[match(overlapping_proteins, top100_list[[m1]]$ID)],
    Target = top100_list[[m1]]$Target[match(overlapping_proteins, top100_list[[m1]]$ID)],
    fdr_m1 = top100_list[[m1]]$fdr[match(overlapping_proteins, top100_list[[m1]]$ID)],
    fdr_m2 = top100_list[[m2]]$fdr[match(overlapping_proteins, top100_list[[m2]]$ID)]
  )
  colnames(result_df)[4:5] <- paste0(method.list, ".pval")
  
  # Prepare Upset plot dataset
  top100_list_clean <- lapply(top100_list, function(x) x$EntrezGeneSymbol[!is.na(x$EntrezGeneSymbol)])
  mlist <- UpSetR::fromList(top100_list_clean) 
  
  # Create UpSet plot
  upset_plot <- upset(mlist, 
                      nsets = length(method.list), 
                      text.scale = c(2, 2, 1.5, 1.5, 2, 2),
                      line.size = 1.5,
                      nintersects = NA, 
                      order.by = "freq")
  
  # Create Venn diagram
  # color_palette <- c("#F7FEF0", "#E9F7E5", "#CEEFCC", "#BFE8C1", "#BCF4C5", "#92C2A6",
  #                    "#D6F6FF", "#ACEEFE", "#6FC8CA", "#58B8D1", "#3492B2", "#04579B")
  color_palette <- c("#F7FEF0", "#BFE8C1", "#58B8D1", "#04579B")
  venn_plot <- venn.diagram(
    x = top100_list_clean,
    filename = NULL,
    fill = color_palette[1:length(method.list)],
    alpha = 0.5,
    cex = 2,
    cat.cex = 2,
    fontfamily = "sans",
    disable.logging = TRUE
  )
  
  return(list(overlapping_proteins = result_df, 
              upset_plot = upset_plot, 
              venn_plot = venn_plot))
}

# Function to perform GO analysis and create plot
perform_go_analysis_and_plot <- function(overlap_final, plotp, method.list) {
  significant_genes <- overlap_final %>% 
    dplyr::select(EntrezGeneSymbol) %>%
    filter(!is.na(EntrezGeneSymbol) & EntrezGeneSymbol != "") %>%
    distinct()
  significant_genes <- bitr(significant_genes$EntrezGeneSymbol, 
                            fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  go_result <- enrichGO(gene = significant_genes$ENTREZID,
                        OrgDb = org.Hs.eg.db,
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.5,
                        qvalueCutoff = 0.5)
  plot <- tryCatch({
    if (nrow(go_result) > 0) {
      dotplot(go_result, showCategory = 15) + 
        ggtitle(paste0("Top ", 20, " GO BP for ", method.list[1]," and ", method.list[2]," Overlap"))
    } else {
      ggplot() + 
        ggtitle("No significant GO Enrichment pathways found") +
        theme_minimal()
    }
  }, error = function(e) {
    message("Error in creating plot: ", e$message)
    ggplot() + 
      ggtitle("Error in creating GO Enrichment pathway plot") +
      theme_minimal()
  })
  return(plot)
}

get_fold_change_matrices <- function(type, brid, prot_cand, prot_interest) {
  if (type == "iA"){
    df_norm <- readRDS("../data/norm_iA.rds")
  } else {
    df_norm <- readRDS("../data/norm_iN.rds")
  }
  
  # adat_file data cleaning
  cleanData <- df_norm %>%
    filter(BRID == brid) %>%
    mutate(bi_group = ifelse(Condition=="WT", 0, 1)) %>% 
    dplyr::select(bi_group, all_of(prot_interest)) %>%
    tidyr::pivot_longer(cols = all_of(prot_interest), 
                        names_to = "ID", values_to = "Expression") %>%
    left_join(prot_cand %>% dplyr::select(Accession, EntrezGeneSymbol, Target),
              by = c("ID" = "Accession")) %>%
    filter(!is.na(EntrezGeneSymbol) & EntrezGeneSymbol != "")
  
  # Calculate fold change for each protein at each timepoint
  fold_change_data <- cleanData %>%
    group_by(ID) %>%
    summarise(
      mean_expr_0 = mean(Expression[bi_group == 0], na.rm = TRUE),
      mean_expr_1 = mean(Expression[bi_group == 1], na.rm = TRUE),
      .groups = "drop") %>%
    mutate(fold_change = log2(mean_expr_1 / mean_expr_0)) %>%
    dplyr::select(ID, fold_change) %>%
    left_join(prot_cand %>% dplyr::select(Accession, EntrezGeneSymbol, Target),
              by = c("ID" = "Accession")) %>%
    arrange(EntrezGeneSymbol)
  
  # Remove columns that contain only NA values
  cols_to_remove <- colnames(fold_change_data)[colSums(is.na(fold_change_data)) == nrow(fold_change_data)]
  fold_change_data <- fold_change_data[, !(colnames(fold_change_data) %in% cols_to_remove)]
  
  return(fold_change_data)
}

# Function to create volcano plot
library(EnhancedVolcano)
volcano_fun <- function(type, brid, volcano_data, padj_thresh = 0.05, FC_thresh = 1.5, max_label = 20) {
  
  # Calculate number of up and down regulated genes
  up_regulated <- sum(volcano_data$padj < padj_thresh & volcano_data$log2FoldChange > log2(FC_thresh), na.rm = TRUE)
  down_regulated <- sum(volcano_data$padj < padj_thresh & volcano_data$log2FoldChange < -log2(FC_thresh), na.rm = TRUE)
  
  # Select top genes to label
  top_genes <- volcano_data %>%
    filter(padj < padj_thresh, abs(log2FoldChange) > log2(FC_thresh)) %>%
    top_n(max_label, wt = abs(log2FoldChange)) %>%
    pull(EntrezGeneSymbol)
  
  p1 <- EnhancedVolcano(volcano_data,
                        lab = volcano_data$EntrezGeneSymbol,
                        x = 'log2FoldChange',
                        y = 'padj',
                        xlim = c(-5, 5),
                        ylim = c(0, 7.5),
                        # title = paste(),
                        subtitle = paste("Volcano Plot -", brid, "KO vs WT.", "Type:", type),
                        caption = paste0("#DEG = ", up_regulated, " up, ", down_regulated, " down"),
                        titleLabSize = 0,
                        subtitleLabSize = 12,
                        captionLabSize = 12,
                        axisLabSize = 12,
                        pCutoff = padj_thresh,
                        FCcutoff = log2(FC_thresh),
                        pointSize = 2,
                        labSize = 2.0,
                        labCol = 'black',
                        labFace = 'bold',
                        boxedLabels = TRUE,
                        colAlpha = 0.4,
                        legendPosition = 'none',
                        # legendLabSize = 10,
                        # legendIconSize = 4.0,
                        drawConnectors = TRUE,
                        widthConnectors = 0.75,
                        colConnectors = 'orange',
                        col=c('black', 'black', 'black', 'red3'),
                        max.overlaps = Inf,
                        selectLab = top_genes)
  
  return(p1)
}
