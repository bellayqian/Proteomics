##The code is comparing btw 
# Condition1_tp1_sample1_r1+r2+r3, with
# Condition2_tp1_sample1_r4+r5+r6.
# .
# So, we are comparing,
# Condition1_tp0_sample1_r1+r2+r3, with
# Condition2_tp0_sample2_r4+r5+r6+r7 +
# Condition2_tp0_sample3_r8+r9+r10+r11 +
# Condition2_tp0_sample4_r12+r13+r14+r15 +
# Condition2_tp0_sample5_r16+r17+r18+r19
# .
# Sample1 contains 1501, 1502, 1503, they are bio replicates
# This is the crude comparison, control group is 1501-1503

getDE<-function(data){
  library(ROTS)
  
  control <- "WT"
  case <- "KO"
  
  temp.data <- data %>% column_to_rownames("protein")
  data1 <- data[, grep(control, colnames(data))]
  data2 <- data[, grep(case, colnames(data))]
  groups=c(rep(1, length(data1)), rep(0, length(data2)))
  
  # Additional check to ensure no rows with all NAs
  all_na_rows <- which(rowSums(!is.na(temp.data)) == 0)
  if (length(all_na_rows) > 0) {
    temp.data <- temp.data[-all_na_rows, ]
  }
    
  rots.out <- ROTS(data = temp.data, groups = groups, B = 100, K = nrow(temp.data) / 5, 
                   paired = F, seed = 1, progress = F, log = T)
  rots.df <- data.frame(p = rots.out$pvalue) #FDR
  rownames(rots.df) <- rownames(rots.out$data)
  
  p.val.df <- rots.df 

  if (!is.data.frame(p.val.df)) {
    p.val.df <- as.data.frame(p.val.df)
  }
  is_single_column <- ncol(p.val.df) == 1
  
  if (is_single_column) {
    p.val.df <- data.frame(
      ProteinID = rownames(p.val.df),
      p.value = p.val.df[,1],
      stringsAsFactors = FALSE
    )
    rownames(p.val.df) <- p.val.df$ProteinID
    rems <- which(is.na(p.val.df$p.value))
    
    p.vals <- p.val.df$p.value
    names(p.vals) <- p.val.df$ProteinID
  }
  
  # Create the result dataframe
  res <- data.frame(
    "Feature ID" = names(p.vals),
    fdr = as.numeric(p.vals),
    stringsAsFactors = FALSE
  )
  
  return(res)
}