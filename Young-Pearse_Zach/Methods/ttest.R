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

getDE<-function(cleanData, tp_ind, des_matrix){

  source("./functions.R")
  control <- "Condition1"
  case <- "Condition2"

  # select_tp <- unique(cleanData$tp_num)[order(unique(cleanData$tp_num))]
  cleanData <- readRDS(paste0(datapath, "_normalized_cleanData.rds"))

  df_processed <- cleanData %>%
    filter(Group %in% c(control_group, treatment_group)) %>%
    mutate(bi_group = ifelse(Group == control_group, 0, 1))

  if (tp_ind != "alltp") {
    print(paste0("ttest: ", tp_ind))
    # Select timepoint
    df_subset <- df_processed %>% filter(tp_num %in% tp_ind)
    table(df_subset$Animal, df_subset$tp_num)

    # Function to perform t-test analysis
    ttest <- function(bi_group, data) {
      t_tests <- getAnalyteInfo(data) %>%
        dplyr::select(AptName, SeqId, Target = TargetFullName, EntrezGeneSymbol, UniProt) %>%
        mutate(
          formula = map(AptName, ~ as.formula(paste(.x, "~ ", "bi_group"))),  # create formula
          t_test  = map(formula, ~ stats::t.test(.x, data = data)),  # fit t-tests
          t_stat  = map_dbl(t_test, "statistic"),                   # extract t-statistic
          p.value = map_dbl(t_test, "p.value"),                     # extract p-values
          fdr     = p.adjust(p.value, method = "BH")                # FDR correction
        ) %>%
        arrange(p.value) %>%       # re-order by p-value
        mutate(rank = row_number()) # add numeric ranks

      # Calculate log2FC
      log2FC <- sapply(t_tests$AptName, function(apt) {
        case_mean <- mean(data[data[[bi_group]] == 1, apt], na.rm = TRUE)
        control_mean <- mean(data[data[[bi_group]] == 0, apt], na.rm = TRUE)
        log2(case_mean / control_mean)
      })

      output <- t_tests %>%
        dplyr::select(AptName, EntrezGeneSymbol, Target, fdr) %>%
        mutate(log2FC = log2FC)

      return(output)
    }

    res <- ttest(bi_group = "bi_group", data = df_subset)
    res <- res[order(res$fdr), ]
  }
  return(res)
}

# For Benchmarking dataset
# getDE <- function(data, tp_ind) {
#   control <- "Condition1"
#   case <- "Condition2"
# 
#   nr.timepoints <- (ncol(data) / 2) / 3 #Applicable for datasets with 2 conditions and 3 replicates in each condition
# 
#   #loop through the timepoints, assumes that the data is ordered according to timepoint and replicates
#   data1 <- data[, grep(control, colnames(data))]
#   data2 <- data[, grep(case, colnames(data))]
# 
#   p.val.df <- data.frame(matrix(nrow = nrow(data1), ncol = nr.timepoints))
#   rownames(p.val.df) <- rownames(data1)
# 
#   #As mentioned above, assumes 3 replicates/condition and the data ordered in a certain way
#   ind1 <- 1
#   ind2 <- 3
# 
#   for (i in 1:nr.timepoints) {
#     dat1 <- data1[, ind1:ind2]
#     dat2 <- data2[, ind1:ind2]
# 
#     temp.data <- cbind(dat1, dat2)
# 
#     #take out proteins with too many missing values for t-test
#     rems1 <- which(rowSums(is.na(dat1)) > 1)
#     rems2 <- which(rowSums(is.na(dat2)) > 1)
#     rems <- c(rems1, rems2)
#     rems <- unique(rems)
#     rem.names <- rownames(temp.data)[rems]
#     if (length(rems) > 0) {
#       temp.data <- temp.data[-rems, ]
#     }
# 
#     ttest.df <- data.frame(p = apply(temp.data, 1, function(x) t.test(x[1:3], x[4:6], na.action = "na.omit")$p.value))
#     rownames(ttest.df) <- rownames(temp.data)
# 
#     if (length(rems) > 0) {
#       rem.frame <- data.frame(matrix(nrow = length(rem.names), ncol = 1))
#       rownames(rem.frame) <- rem.names
#       colnames(rem.frame) <- "p"
#       ttest.df <- rbind(ttest.df, rem.frame)
#     }
# 
#     p.val.df[, i] <- ttest.df[match(rownames(p.val.df), rownames(ttest.df)), ]
# 
#     ind1 <- ind1 + 3
#     ind2 <- ind2 + 3
#   }
# 
#   rems <- which(apply(p.val.df, 1, function(x) length(which(is.na(x)))) == nr.timepoints)
#   rem.names <- rownames(p.val.df)[rems]
#   rems2 <- rep(NA, length(rems))
#   names(rems2) <- rem.names
# 
#   if (length(rems) > 0) {
#     p.val.df <- p.val.df[-rems, ]
#   }
# 
#   p.vals <- apply(p.val.df, 1, function(x) min(x, na.rm = T)) #Summarize with minimum
#   names(p.vals) <- rownames(p.val.df)
# 
#   if (length(rems) > 0) {
#     p.vals <- c(p.vals, rems2)
#   }
# 
#   res <- cbind(names(p.vals), p.vals)
#   res[is.nan(res)] <- NA
# 
#   res <- data.frame(res, stringsAsFactors = FALSE)
#   res[, 1] <- as.character(res[, 1])
#   res[, 2] <- as.numeric(res[, 2])
#   colnames(res) <- c("ID", "pval")
# 
#   if (any(is.nan(res[, 2]))) {
#     res[which(is.nan(res[, 2])), 2] <- NA
#   }
# 
#   return(res)
# }