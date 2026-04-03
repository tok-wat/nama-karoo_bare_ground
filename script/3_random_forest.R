library(tidyverse)
library(collinear)
library(caret)
library(UBL)

##################################
rm(list=ls())

gc(reset = TRUE)
gc(reset = TRUE)
##################################
# Load a script
source("./Script/function/functions_util.R")
source("./Script/function/functions_rf.R")
##################################

# ------------------------------------------------------------------------------
# Load training/validation split df files and specifications of predictor variables
sf_growing_split <- readRDS("./Script/tmp/sf_growing_split.rds")
sf_non_growing_split <- readRDS("./Script/tmp/sf_non_growing_split.rds")

# ------------------------------------------------------------------------------
# Running RF pipelines with different datasets and preprocessing methods

rf_growing <- func_rf_pipeline(dataframe = sf_growing_split$train_df,
                 predictors = c(vi),
                 blocks = sf_growing_split$spatial_cv_block,
                 seed=1234)

rf_growing_env <- func_rf_pipeline(dataframe = sf_growing_split$train_df,
                 predictors = c(vi, soil,minerals,topography),
                 blocks = sf_growing_split$spatial_cv_block,
                 seed=1234)

rf_growing_smote <- func_rf_smote_pipeline(dataframe = sf_growing_split$train_df,
                       predictors = c(vi),
                       blocks = sf_growing_split,
                       seed=1234)

rf_growing_env_smote <- func_rf_smote_pipeline(dataframe = sf_growing_split$train_df,
                       predictors = c(vi, soil,minerals,topography),
                       blocks = sf_growing_split,
                       seed=1234)

rf_non_growing <- func_rf_pipeline(dataframe = sf_non_growing_split$train_df,
                 predictors = c(vi),
                 blocks = sf_non_growing_split$spatial_cv_block,
                 seed=1234)

rf_non_growing_env <-func_rf_pipeline(dataframe = sf_non_growing_split$train_df,
                 predictors = c(vi, soil,minerals,topography),
                 blocks = sf_non_growing_split$spatial_cv_block,
                 seed=1234)

# ------------------------------------------------------------------------------
# Save trained RF models to .rds files

# Growing season models
saveRDS(rf_growing$best_model,           file = "./Script/models/growing/rf_1-1_rf_growing.rds")
saveRDS(rf_growing_env$best_model,       file = "./Script/models/growing/rf_1-2_rf_growing_env.rds")
saveRDS(rf_growing_smote$best_model,     file = "./Script/models/growing/rf_1-3_rf_growing_smote.rds")
saveRDS(rf_growing_env_smote$best_model, file = "./Script/models/growing/rf_1-4_rf_growing_env_smote.rds")

# Non-growing season models
saveRDS(rf_non_growing$best_model,       file = "./Script/models/non_growing/rf_2-1_rf_non_growing.rds")
saveRDS(rf_non_growing_env$best_model,   file = "./Script/models/non_growing/rf_2-2_rf_non_growing_env.rds")

# ------------------------------------------------------------------------------
# Functions to extract model summary (best performance + selected predictors)

# For models without SMOTE
rf_summary <- function(rf_model) {
  rf_model$best_model$results %>%
    filter(mtry == rf_model$best_model$bestTune$mtry) %>%
    mutate(predictors = paste0(rf_model$RFE_output$selected_predictors, collapse = ", "))
}

# For models with SMOTE
rf_smote_summary <- function(rf_smote_model) {
  rf_smote_model$cv_output %>%
    mutate(mtry = rf_smote_model$best_model$bestTune$mtry) %>%
    mutate(predictors = paste0(rf_smote_model$RFE_output$selected_predictors, collapse = ", "))
}

# ------------------------------------------------------------------------------
# Combine summaries from all models

summary_rf_model <- rbind(
  rf_summary(rf_growing),
  rf_summary(rf_growing_env),
  rf_smote_summary(rf_growing_smote),
  rf_smote_summary(rf_growing_env_smote),
  rf_summary(rf_non_growing),
  rf_summary(rf_non_growing_env)
) %>%
  mutate(
    model = paste0("rf_", c(rep(1, 4), rep(2, 2)), "-", c(1:4, 1:2)),
    season = c(rep("growing", 4), rep("non-growing", 2)),
    preprocess = c("none", "none", "smote", "smote", "none", "none")
  )

# Convert to character matrix before saving (if needed)
summary_rf_model <- apply(summary_rf_model, 2, as.character)

# ------------------------------------------------------------------------------
# Save model summary to CSV
write.csv(summary_rf_model, "./Script/tmp/rf_summary.csv", row.names = FALSE)
