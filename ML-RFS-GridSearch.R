####################################################
# RFS WITH HYPERPARAMETER TUNING 
####################################################

# Load required libraries without startup messages
suppressPackageStartupMessages({
  library(survival)
  library(readxl)
  library(dplyr)
  library(riskRegression)
  library(ggplot2)
  library(caret)      # For cross-validation folds
  library(pec)        # For cindex and calPlot functions
  library(randomForestSRC)
  library(tidyr)      # For data manipulation
})

#--------------------------------------------------------------------------
# Step 1: Load the Data
#--------------------------------------------------------------------------
file_path <- "path_to_your_file"
df <- read_excel(file_path, sheet = "Factors") %>%
  as.data.frame()

#--------------------------------------------------------------------------
# Step 2: Prepare the Data
#--------------------------------------------------------------------------
# Convert categorical variables to factors
categorical_vars <- c('Sex', 'ECOG_Merged', 'Smoking_Type', 'Drinking_History_Merged',
                      'T', 'Node', 'Primary_Treatment_Modality', 'Insurance_Type',
                      'Periodontal_Grading', 'HPV', 'Extraction_befor_RT_End',
                      'Disease_Site_Merged_2', 'Histological_Diagnosis', 'Chemotherapy')
df[categorical_vars] <- lapply(df[categorical_vars], as.factor)

# Define continuous variables
continuous_vars <- c('Age', 'Smoking_Pack_per_Year', 'DMFS', 'Income_1000',
                     'Number_Extracted_Teeth', 'Number_Teeth_after_Extraction',
                     'Duration_RT_Days', 'RT_Dose', 'D20')

# Ensure survival-related variables are numeric
df$time <- as.numeric(as.character(df$ClinRad_Time_Indicator_M...8))
df$delta <- as.numeric(as.character(df$ClinRad_M_Competing))  # 0: censored, 1: event, 2: competing risk

# Remove rows with missing values in predictors or outcome
df <- df %>% 
  select(all_of(c(categorical_vars, continuous_vars, "time", "delta"))) %>% 
  na.omit()

# Cap time at 114 months
df <- df %>% mutate(time = ifelse(time > 114, 114, time))

#--------------------------------------------------------------------------
# Step 3: Set up Repeated 5-Fold Cross-Validation (5×5 = 25 splits)
#--------------------------------------------------------------------------
set.seed(123)  # For reproducibility
k_folds <- 5
n_repeats <- 5

# createMultiFolds() returns a named list of training indices for each fold × repeat
multi_folds <- createMultiFolds(factor(df$delta), k = k_folds, times = n_repeats)
all_indices <- seq_len(nrow(df))

#--------------------------------------------------------------------------
# Step 4: Define Hyperparameter Grid
#--------------------------------------------------------------------------
# Only tuning nodesize and nodedepth
hyper_grid <- expand.grid(
  nodesize = c(15),    # Example values; adjust as needed
  nodedepth = c(20)   # Example values; adjust as needed
)

#--------------------------------------------------------------------------
# Step 5: Initialize Data Structures to Store Results
#--------------------------------------------------------------------------
# To store performance metrics for each hyperparameter combination
performance_results <- data.frame()

#--------------------------------------------------------------------------
# Step 6: Define Helper Function for Computing Mean and 95% CI
#--------------------------------------------------------------------------
compute_mean_ci <- function(x, ci_level = 0.95) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 2) {
    return(c(Mean = mean(x), Lower = NA, Upper = NA))
  }
  m <- mean(x)
  sd_val <- sd(x)
  se <- sd_val / sqrt(n)
  alpha <- 1 - ci_level
  z <- qnorm(1 - alpha/2)  # ~1.96 for 95% CI
  lower <- m - z * se
  upper <- m + z * se
  c(Mean = m, Lower = lower, Upper = upper)
}

# Specify time horizons of interest
time_horizons <- c(60, 114)

#--------------------------------------------------------------------------
# Step 7: Hyperparameter Tuning Loop
#--------------------------------------------------------------------------
# Iterate over each hyperparameter combination
for (i in 1:nrow(hyper_grid)) {
  current_hyper <- hyper_grid[i, ]
  cat("\nEvaluating Hyperparameters:", 
      "nodesize =", current_hyper$nodesize, 
      ", nodedepth =", current_hyper$nodedepth, "\n")
  
  # Initialize vectors to store IBS and C-index for all folds
  ibs_values <- c()
  cindex_values <- c()
  
  # Iterate over each of the 25 folds
  for (fold_name in names(multi_folds)) {
    cat("  Processing Fold:", fold_name, "\n")
    
    # Get training and test indices
    train_indices <- multi_folds[[fold_name]]
    test_indices  <- setdiff(all_indices, train_indices)
    
    train_data <- df[train_indices, ]
    test_data  <- df[test_indices, ]
    
    # Build the formula for the Fine-Gray model
    predictor_vars <- paste(c(categorical_vars, continuous_vars), collapse = " + ")
    rfs_formula <- as.formula(paste0("Surv(time, delta) ~ ", predictor_vars))
    
    # Fit the Random Forest model with current hyperparameters
    rfs_model <- rfsrc(
      formula = rfs_formula,
      data = train_data,
      ntree = 300,                # Keeping ntree fixed; adjust if desired
      nodesize = current_hyper$nodesize,
      nodedepth = current_hyper$nodedepth,
      splitrule = 'logrankCR',
      cause = 1
    )
    
    # Predict CIF on the test set
    time_points <- seq(0, 114, by = 1)
    test_cif <- predictRisk(rfs_model, newdata = test_data, times = time_points, cause = 1)
    
    # Evaluate performance metrics at specified time horizons
    for (time_horizon in time_horizons) {
      cat("    Evaluating up to time horizon:", time_horizon, "\n")
      
      time_seq <- seq(1, time_horizon, by = 1)
      
      scores <- Score(
        object       = list("RFS" = rfs_model),  
        formula      = Hist(time, delta) ~ 1,
        data         = test_data,
        cause        = 1,
        summary      = "ibs",  # to compute IBS
        times        = time_seq,
        metrics      = c("AUC", "Brier"),
        cens.model   = "km",
        split.method = "none",
        verbose      = FALSE
      )
      
      # Extract IBS
      ibs_value <- scores$Brier$score %>%
        filter(model == "RFS", times == max(time_seq)) %>%
        pull(IBS)
      if (length(ibs_value) == 0) ibs_value <- NA
      ibs_values <- c(ibs_values, ibs_value)
      
      # Extract C-index
      cindex_result <- pec::cindex(
        object      = list("RFS" = rfs_model),
        formula     = Hist(time, delta) ~ 1,
        data        = test_data,
        eval.times  = time_horizon,
        cause       = 1,
        cens.model  = "marginal",
        splitMethod = "none",
        keep.matrix = FALSE,
        verbose     = FALSE
      )
      cindex_value <- cindex_result$AppCindex$RFS
      cindex_values <- c(cindex_values, cindex_value)
    }
  }
  
  # Compute average IBS and C-index across all folds and time horizons
  avg_ibs <- mean(ibs_values, na.rm = TRUE)
  avg_cindex <- mean(cindex_values, na.rm = TRUE)
  
  # Store the results
  performance_results <- rbind(performance_results, 
                                data.frame(
                                  nodesize = current_hyper$nodesize,
                                  nodedepth = current_hyper$nodedepth,
                                  Average_IBS = avg_ibs,
                                  Average_Cindex = avg_cindex
                                ))
}

#--------------------------------------------------------------------------
# Step 8: Identify the Best Hyperparameter Combination
#--------------------------------------------------------------------------
# Define the criterion for "best". For example, lowest IBS and highest C-index.
# Here, we'll prioritize the lowest Average_IBS. In case of ties, the highest Cindex is chosen.
best_hyper <- performance_results %>%
  arrange(Average_IBS, desc(Average_Cindex)) %>%
  slice(1)

cat("\nBest Hyperparameters based on Average IBS and C-index:\n")
print(best_hyper)

#--------------------------------------------------------------------------
# Step 9: Retrain Models with Best Hyperparameters and Collect Metrics
#--------------------------------------------------------------------------
# Initialize lists to store performance metrics
auc_df_list      <- list()
brier_df_list    <- list()
ibs_values_list  <- list()
cindex_df_list   <- list()
cif_list         <- list()

# Extract best hyperparameters
best_nodesize <- best_hyper$nodesize
best_nodedepth <- best_hyper$nodedepth
cat("\nRetraining models with Best Hyperparameters:",
    "nodesize =", best_nodesize, 
    ", nodedepth =", best_nodedepth, "\n")

# Iterate over each of the 25 folds with best hyperparameters
for (fold_name in names(multi_folds)) {
  cat("Processing Fold:", fold_name, "\n")
  
  # Get training and test indices
  train_indices <- multi_folds[[fold_name]]
  test_indices  <- setdiff(all_indices, train_indices)
  
  train_data <- df[train_indices, ]
  test_data  <- df[test_indices, ]
  
  # Build the formula for the Fine-Gray model
  predictor_vars <- paste(c(categorical_vars, continuous_vars), collapse = " + ")
  rfs_formula <- as.formula(paste0("Surv(time, delta) ~ ", predictor_vars))
  
  # Fit the Random Forest model with best hyperparameters
  rfs_model <- rfsrc(
    formula = rfs_formula,
    data = train_data,
    ntree = 300,                # Keeping ntree fixed; adjust if desired
    nodesize = best_nodesize,
    nodedepth = best_nodedepth,
    splitrule = 'logrankCR',
    cause = 1
  )
  
  # Predict CIF on the test set
  time_points <- seq(0, 114, by = 1)
  test_cif <- predictRisk(rfs_model, newdata = test_data, times = time_points, cause = 1)
  
  # Average CIF across test patients for each time point
  cif_df <- data.frame(
    Time = time_points,
    CIF  = colMeans(test_cif),
    Fold = fold_name
  )
  cif_list[[fold_name]] <- cif_df
  
  # Evaluate performance metrics at specified time horizons
  for (time_horizon in time_horizons) {
    cat("  Evaluating up to time horizon:", time_horizon, "\n")
    
    time_seq <- seq(1, time_horizon, by = 1)
    
    scores <- Score(
      object       = list("RFS" = rfs_model),  
      formula      = Hist(time, delta) ~ 1,
      data         = test_data,
      cause        = 1,
      summary      = "ibs",  # to compute IBS
      times        = time_seq,
      metrics      = c("AUC", "Brier"),
      cens.model   = "km",
      split.method = "none",
      verbose      = FALSE
    )
    
    # AUC
    auc_values <- scores$AUC$score %>%
      filter(model == "RFS") %>%
      select(times, AUC)
    auc_values$Fold <- fold_name
    auc_values$Time_Horizon <- time_horizon
    
    # Brier
    brier_values <- scores$Brier$score %>%
      filter(model == "RFS") %>%
      select(times, Brier)
    brier_values$Fold <- fold_name
    brier_values$Time_Horizon <- time_horizon
    
    # IBS for the last time in time_seq
    ibs_value <- scores$Brier$score %>%
      filter(model == "RFS", times == max(time_seq)) %>%
      pull(IBS)
    if (length(ibs_value) == 0) ibs_value <- NA
    
    ibs_value_df <- data.frame(
      Fold         = fold_name,
      Time_Horizon = time_horizon,
      IBS          = ibs_value
    )
    
    # C-index
    cindex_result <- pec::cindex(
      object      = list("RFS" = rfs_model),
      formula     = Hist(time, delta) ~ 1,
      data        = test_data,
      eval.times  = time_horizon,
      cause       = 1,
      cens.model  = "marginal",
      splitMethod = "none",
      keep.matrix = FALSE,
      verbose     = FALSE
    )
    cindex_value <- cindex_result$AppCindex$RFS
    cindex_df <- data.frame(
      Fold         = fold_name,
      Time_Horizon = time_horizon,
      Cindex       = cindex_value
    )
    
    # Store metrics
    auc_df_list[[length(auc_df_list) + 1]]     <- auc_values
    brier_df_list[[length(brier_df_list) + 1]] <- brier_values
    ibs_values_list[[length(ibs_values_list) + 1]] <- ibs_value_df
    cindex_df_list[[length(cindex_df_list) + 1]] <- cindex_df
  }
}

#--------------------------------------------------------------------------
# Step 10: Combine Results (All 25 Splits) into Data Frames
#--------------------------------------------------------------------------
auc_df    <- do.call(rbind, auc_df_list)
brier_df  <- do.call(rbind, brier_df_list)
ibs_df    <- do.call(rbind, ibs_values_list)
cindex_df <- do.call(rbind, cindex_df_list)
cif_all   <- do.call(rbind, cif_list)

#--------------------------------------------------------------------------
# Step 11: Identify the Best Model Based on C-index at Time Horizon 114
#--------------------------------------------------------------------------
cindex_114 <- cindex_df %>%
  filter(Time_Horizon == 114)

best_fold <- cindex_114 %>%
  arrange(desc(Cindex)) %>%
  slice(1) %>%
  pull(Fold)

cat("\nBest fold (model) is Fold:", best_fold, 
    "with C-index at 114 =", 
    cindex_114 %>% filter(Fold == best_fold) %>% pull(Cindex), "\n")

#--------------------------------------------------------------------------
# Step 12: Compute Summary Statistics (Mean + 95% CI) for Metrics
#--------------------------------------------------------------------------
## AUC summary
auc_summary <- auc_df %>%
  group_by(times, Time_Horizon) %>%
  summarise(
    AUC_Mean    = compute_mean_ci(AUC)[1],
    AUC_LowerCI = compute_mean_ci(AUC)[2],
    AUC_UpperCI = compute_mean_ci(AUC)[3],
    .groups     = "drop"
  )

## Brier summary
brier_summary <- brier_df %>%
  group_by(times, Time_Horizon) %>%
  summarise(
    Brier_Mean    = compute_mean_ci(Brier)[1],
    Brier_LowerCI = compute_mean_ci(Brier)[2],
    Brier_UpperCI = compute_mean_ci(Brier)[3],
    .groups       = "drop"
  )

## IBS summary
ibs_summary <- ibs_df %>%
  group_by(Time_Horizon) %>%
  summarise(
    IBS_Mean    = compute_mean_ci(IBS)[1],
    IBS_LowerCI = compute_mean_ci(IBS)[2],
    IBS_UpperCI = compute_mean_ci(IBS)[3],
    .groups     = "drop"
  )

## C-index summary
cindex_summary <- cindex_df %>%
  group_by(Time_Horizon) %>%
  summarise(
    Cindex_Mean    = compute_mean_ci(Cindex)[1],
    Cindex_LowerCI = compute_mean_ci(Cindex)[2],
    Cindex_UpperCI = compute_mean_ci(Cindex)[3],
    .groups        = "drop"
  )

# Categorize time ranges in AUC and Brier data frames (for plotting)
auc_summary <- auc_summary %>%
  mutate(Time_Range = case_when(
    times > 0 & times <= 60 ~ "0-60",
    times > 60 ~ "60-Max"
  ))

brier_summary <- brier_summary %>%
  mutate(Time_Range = case_when(
    times > 0 & times <= 60 ~ "0-60",
    times > 60 ~ "60-Max"
  ))

#--------------------------------------------------------------------------
# Step 13: Plotting with 95% CI Ribbons
#--------------------------------------------------------------------------
# A) Time-dependent AUC
ggplot(auc_summary, aes(x = times, y = AUC_Mean, color = Time_Range)) +
  geom_line() +
  geom_point(size = 0.5) +
  geom_ribbon(aes(ymin = AUC_LowerCI, ymax = AUC_UpperCI, fill = Time_Range),
              alpha = 0.2, color = NA) +
  ggtitle("Repeated 5×5 CV Time-dependent AUC with Time Ranges (95% CI)") +
  xlab("Time (months)") +
  ylab("AUC") +
  scale_x_continuous(
    breaks = seq(0, max(auc_summary$times, na.rm = TRUE), by = 10),
    limits = c(0, max(auc_summary$times, na.rm = TRUE))
  ) +
  scale_color_manual(
    values = c("0-60" = "green",  "60-Max" = "red"),
    name = "Time Range"
  ) +
  scale_fill_manual(
    values = c("0-60" = "green",  "60-Max" = "red"),
    name = "Time Range"
  ) +
  theme_minimal()

# B) Time-dependent Brier score
ggplot(brier_summary, aes(x = times, y = Brier_Mean, color = Time_Range)) +
  geom_line() +
  geom_point(size = 0.5) +
  geom_ribbon(aes(ymin = Brier_LowerCI, ymax = Brier_UpperCI, fill = Time_Range),
              alpha = 0.2, color = NA) +
  ggtitle("Repeated 5×5 CV Brier Scores with Time Ranges (95% CI)") +
  xlab("Time (months)") +
  ylab("Brier Score") +
  scale_x_continuous(
    breaks = seq(0, max(brier_summary$times, na.rm = TRUE), by = 10),
    limits = c(0, max(brier_summary$times, na.rm = TRUE))
  ) +
  scale_color_manual(
    values = c("0-60" = "green",  "60-Max" = "red"),
    name = "Time Range"
  ) +
  scale_fill_manual(
    values = c("0-60" = "green",  "60-Max" = "red"),
    name = "Time Range"
  ) +
  theme_minimal()

#--------------------------------------------------------------------------
# Step 14: Print IBS and C-index Summaries
#--------------------------------------------------------------------------
cat("\nRepeated (5×5) Cross-Validated Integrated Brier Scores (IBS) for each time horizon:\n")
print(ibs_summary)

cat("\nRepeated (5×5) Cross-Validated C-index for each time horizon:\n")
print(cindex_summary)

#--------------------------------------------------------------------------
# Step 15: Combine & Summarize CIF
#--------------------------------------------------------------------------
cif_summary <- cif_all %>%
  group_by(Time) %>%
  summarise(
    CIF_Mean    = compute_mean_ci(CIF)[1],
    CIF_LowerCI = compute_mean_ci(CIF)[2],
    CIF_UpperCI = compute_mean_ci(CIF)[3],
    .groups     = "drop"
  )

# Extract CIF values for months 60 and 114 (if present)
cif_at_60  <- cif_summary %>% filter(Time == 60)  %>% pull(CIF_Mean)
cif_at_114 <- cif_summary %>% filter(Time == 114) %>% pull(CIF_Mean)
cat("\nAverage CIF at month 60:", cif_at_60, "\n")
cat("Average CIF at month 114:", cif_at_114, "\n")

#--------------------------------------------------------------------------
# Step 16: Plot the Average CIF with 95% CI
#--------------------------------------------------------------------------
ggplot(cif_summary, aes(x = Time, y = CIF_Mean)) +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = CIF_LowerCI, ymax = CIF_UpperCI),
              fill = "blue", alpha = 0.2) +
  geom_point(size = 1, color = "blue") +
  ggtitle("Repeated (5×5) CV Average CIF Across Folds (95% CI)") +
  xlab("Time (months)") +
  ylab("Average CIF") +
  scale_x_continuous(
    breaks = seq(0, 114, by = 10),
    limits = c(0, 114)
  ) +
  theme_minimal()

#--------------------------------------------------------------------------
# Step 17: Summarize Performance Results for Hyperparameter Combinations
#--------------------------------------------------------------------------
# Compute the average C-index and IBS for each hyperparameter combination
summary_table <- performance_results %>%
  group_by(nodesize, nodedepth) %>%
  summarise(
    Mean_IBS = mean(Average_IBS, na.rm = TRUE),
    Mean_Cindex = mean(Average_Cindex, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Mean_IBS, desc(Mean_Cindex))

# Print the summary table
cat("\nSummary of Average IBS and C-index for each Hyperparameter Combination:\n")
print(summary_table)


```
