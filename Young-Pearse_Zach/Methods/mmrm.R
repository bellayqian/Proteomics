# install.packages("mmrm")
# load("UPS1_Mix_Full_ExampleData.RData")
# load("UPSSpikeProteins.RData")

getDE<-function(data, tp_ind, des_matrix){

  library(mmrm)
  library(emmeans)
  library(dplyr)
  library(tidyr)
  library(data.table)
  source("./functions.R")
  control <- "Condition1"
  case <- "Condition2"

  data <- as.data.frame(data) %>% rownames_to_column(var = "Protein")
  
  if (oligoID == "oligo368"){
    long_data <- data %>%
      pivot_longer(cols = -Protein,
                   names_to = c("Condition", "Timepoint", "Sample", "Replicate"),
                   names_pattern = "(Condition\\d)_timepoint_(\\d+)_Sample_(\\d+)_r(\\d)",
                   values_to = "Expression") %>%
      mutate(Treatment = case_when(Condition == "Condition1" ~ "Control",
                                   Condition == "Condition2" ~ "Treatment"),
             Timepoint = as.numeric(Timepoint)) %>%
      group_by(Protein, Sample) %>%
      mutate(responseBL = Expression[which(Timepoint == 0)[1]],
             response = if_else(Timepoint != 0, Expression, NA_real_)) %>%
      ungroup() %>%
      filter(!(Timepoint %in% c(0, 1))) %>%
      dplyr::select(Protein, subject = Sample, trt = Treatment,
                    time = Timepoint, response, responseBL) %>%
      mutate(subject = factor(subject),
             trt = factor(trt),
             time = factor(time)) %>%
      arrange(subject, Protein)
  } else {
    long_data <- data %>%
      pivot_longer(cols = -Protein,
                   names_to = c("Condition", "Timepoint", "Sample", "Replicate"),
                   names_pattern = "(Condition\\d)_timepoint_(\\d+)_Sample_(\\d+)_r(\\d)",
                   values_to = "Expression") %>%
      mutate(Treatment = case_when(Condition == "Condition1" ~ "Control",
                                   Condition == "Condition2" ~ "Treatment"),
             Timepoint = as.numeric(Timepoint)) %>%
      group_by(Protein, Sample) %>%
      mutate(responseBL = Expression[which(Timepoint == 0)[1]],
             response = if_else(Timepoint != 0, Expression, NA_real_)) %>%
      ungroup() %>%
      filter(!(Timepoint %in% c(0))) %>%
      dplyr::select(Protein, subject = Sample, trt = Treatment,
                    time = Timepoint, response, responseBL) %>%
      mutate(subject = factor(subject),
             trt = factor(trt),
             time = factor(time)) %>%
      arrange(subject, Protein)
  }

  # set.seed(1)
  # selected_proteins <- sample(unique(long_data$Protein), 100)

  # Filter the long_data to include only the selected proteins
  # long_data_sample <- long_data %>% filter(Protein %in% selected_proteins)

  # Get unique proteins
  proteins <- unique(long_data$Protein)

  results_df <- data.frame(Protein = character(), pval = numeric(), stringsAsFactors = FALSE)

  # Loop through each protein
  for (protein_name in proteins) {

    # Subset data for the current protein
    protein_data <- long_data %>%
      filter(Protein == protein_name)

    # Fit the MMRM model
    fit.model <- tryCatch({
      mmrm(formula = response ~ responseBL + trt + time + trt:time + us(time|subject),
           data = protein_data,
           method = "Kenward-Roger",
           vcov = "Kenward-Roger-Linear")
    }, error = function(e) {
      tryCatch({
        mmrm(formula = response ~ responseBL + trt + time + trt:time + toeph(time|subject),
             data = protein_data,
             method = "Kenward-Roger",
             vcov = "Kenward-Roger-Linear")
      }, error = function(e) {
        # message(paste("Error in fitting both models for protein:", protein_name))
        return(NULL)
      })
    })

    if (!is.null(fit.model)) {
      # Calculate least squares means
      # emm <- emmeans(fit.model, specs = c("trt", "time"))

      # Within treatment differences
      # within_trt_diffs <- pairs(lsmeans, by = "trt", reverse = TRUE)
      # confint(within_trt_diffs)

      # Between treatment differences
      # between_trt_diffs  <- pairs(lsmeans, by = "time", reverse = TRUE)
      # confint(between_trt_diffs)

      # Longitudinal analysis
      # Perform ANOVA to get overall treatment effect
      anova_result <- tryCatch({
        car::Anova(fit.model)
      }, error = function(e) {

        # Try fitting the model with 'toeph' covariance structure
        fit.model_toeph <- tryCatch({
          mmrm(formula = response ~ responseBL + trt + time + trt:time + toeph(time|subject),
               data = protein_data,
               method = "Kenward-Roger",
               vcov = "Kenward-Roger-Linear")
        }, error = function(e2) {
          return(NULL)
        })

        if (!is.null(fit.model_toeph)) {
          tryCatch({
            Anova(fit.model_toeph)
          }, error = function(e3) {
            return(NULL)
          })
        } else {
          return(NULL)
        }
      })

      if (!is.null(anova_result)) {
        # Extract p-value for the treatment effect (main effect of trt)
        p_value <- anova_result["trt", "Pr(>=F)"]

        new_row <- data.frame(Protein = protein_name, pval = p_value)
        results_df <- rbind(results_df, new_row)
      }
    }
  }

  res <- results_df %>%
    # mutate(adj_pval = p.adjust(pval, method = "fdr")) %>%
    dplyr::select("ID" = Protein, "P value" = pval) %>%
    arrange(`P value`)

  return(res)
}

# # Benchmarking
# getDE<-function(data, tp_ind){
#   
#   library(mmrm)
#   library(emmeans)
#   library(car)
#   library(data.table)
#   source("./functions.R")
#   control <- "Condition1"
#   case <- "Condition2"
#   
#   load("UPS1_Mix_Full_ExampleData.RData")
#   load("UPSSpikeProteins.RData")
#   
#   data <- as.data.frame(data) %>% rownames_to_column(var = "Protein")
#   
#   long_data <- data %>%
#     pivot_longer(cols = -Protein,
#                  names_to = "Sample",
#                  values_to = "Expression") %>%
#     mutate(Condition = str_extract(Sample, "Condition\\d"),
#            Treatment = case_when(Condition == "Condition1" ~ "Control",
#                                  Condition == "Condition2" ~ "Treatment"),
#            Timepoint = as.numeric(str_extract(Sample, "(?<=timepoint_)\\d+"))) %>%
#     # group_by(Protein, Sample) %>%
#     # mutate(responseBL = Expression[which(Timepoint == 1)[1]]) %>%
#     # ungroup() %>%
#     # filter(!(Timepoint %in% c(1))) %>%
#     dplyr::select(Protein, subject = Sample, trt = Treatment, 
#                   time = Timepoint, response = Expression) %>%
#     mutate(subject = factor(subject),
#            trt = factor(trt),
#            time = factor(time)) %>%
#     arrange(subject, Protein)
#   
#   # set.seed(1)
#   # selected_proteins <- sample(unique(long_data$Protein), 100)
#   # # Filter the long_data to include only the selected proteins
#   # long_data_sample <- long_data %>% filter(Protein %in% selected_proteins)
#   
#   # Get unique proteins
#   proteins <- unique(long_data$Protein)
#   
#   results_df <- data.frame(Protein = character(), pval = numeric(), stringsAsFactors = FALSE)
#   
#   # Loop through each protein
#   for (protein_name in proteins) {
#     
#     # Subset data for the current protein
#     protein_data <- long_data %>%
#       filter(Protein == protein_name) #"O00762ups|UBE2C_HUMAN_UPS")
#     
#     # Fit the MMRM model
#     fit.model <- tryCatch({
#       mmrm(formula = response ~ trt + time + trt:time + us(time|subject),
#            data = protein_data,
#            method = "Kenward-Roger",
#            vcov = "Kenward-Roger-Linear")
#     }, error = function(e) {
#       tryCatch({
#         mmrm(formula = response ~ trt + time + trt:time + toeph(time|subject),
#              data = protein_data,
#              method = "Kenward-Roger",
#              vcov = "Kenward-Roger-Linear")
#       }, error = function(e) {
#         message(paste("Error in fitting both models for protein:", protein_name))
#         return(NULL)
#       })
#     })
#     
#     if (!is.null(fit.model)) {
#       # Calculate least squares means
#       # emm <- emmeans(fit.model, specs = c("trt", "time"))
#       
#       # Within treatment differences
#       # within_trt_diffs <- pairs(lsmeans, by = "trt", reverse = TRUE)
#       # confint(within_trt_diffs)
#       
#       # Between treatment differences
#       # between_trt_diffs  <- pairs(lsmeans, by = "time", reverse = TRUE)
#       # confint(between_trt_diffs)
#       
#       # Longitudinal analysis
#       # Perform ANOVA to get overall treatment effect
#       anova_result <- tryCatch({
#         car::Anova(fit.model)
#       }, error = function(e) {
#         
#         # Try fitting the model with 'toeph' covariance structure
#         fit.model_toeph <- tryCatch({
#           mmrm(formula = response ~ trt + time + trt:time + toeph(time|subject),
#                data = protein_data,
#                method = "Kenward-Roger",
#                vcov = "Kenward-Roger-Linear")
#         }, error = function(e2) {
#           return(NULL)
#         })
#         
#         if (!is.null(fit.model_toeph)) {
#           tryCatch({
#             Anova(fit.model_toeph)
#           }, error = function(e3) {
#             return(NULL)
#           })
#         } else {
#           return(NULL)
#         }
#       })
#       
#       if (!is.null(anova_result)) {
#         # Extract p-value for the treatment effect (main effect of trt)
#         p_value <- anova_result["trt", "Pr(>=F)"]
#         
#         new_row <- data.frame(Protein = protein_name, pval = p_value)
#         results_df <- rbind(results_df, new_row)
#       }
#     }
#   }
#   
#   res <- results_df %>% 
#     # mutate(adj_pval = p.adjust(pval, method = "fdr")) %>%
#     dplyr::select("ID" = Protein, "P value" = pval) %>%
#     arrange(`P value`)
#   
#   return(res)
# }