##Assumes that the samples in the given data are organized as:
#Condition 1 timepoint 1 replicate 1,
#Condition 1 timepoint 1 replicate 2,
#Condition 1 timepoint 1 replicate 3,
#Condition 1 timepoint 2 replicate 1,
#Condition 1 timepoint 2 replicate 2,
#.
#.
#.
#Condition 2 timepoint 1 replicate 1,
#Condition 2 timepoint 1 replicate 2,
#Condition 2 timepoint 1 replicate 3,
#Condition 2 timepoint 2 replicate 1,
#Condition 2 timepoint 2 replicate 2,

getDE<-function(data){
  library(Biobase)
  library(limma)
  
  control <- "WT"
  case <- "KO"
  
  data <- data %>% column_to_rownames("protein")
  
  data1 <- data[, grep(control, colnames(data))]
  data2 <- data[, grep(case, colnames(data))]
  data.temp <- data[, c(grep(case, colnames(data)), grep(control, colnames(data)))]
  
  case.group <- rep("BR24_KO", ncol(data2))
  control.group <- rep("BR24_WT", ncol(data1))
  
  f <- factor(c(case.group, control.group))
  design <- model.matrix(~0 + f)
  colnames(design) <- levels(f)
  fit <- lmFit(data.temp, design)
  
  mc <- makeContrasts(contrasts = "BR24_KO - BR24_WT", levels = design)
  c.fit <- contrasts.fit(fit, mc)
  eb <- eBayes(c.fit)
  
  top_proteins <- topTable(eb, number = Inf, sort.by = "P")
  
  res <- top_proteins %>% 
    dplyr::select(logFC, adj.P.Val) %>% 
    rownames_to_column("Feature.ID") %>% 
    arrange(logFC, adj.P.Val)

  return(res)
}
