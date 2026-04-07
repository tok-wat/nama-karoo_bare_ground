# Load required libraries

# data wrangling
library(tidyverse)
library(jsonlite)

# data analysis
library(sf)
library(MASS)  # for box-cox transformation
library(CAST)  # for knndm and CreateSpacetimeFolds
library(caret) # wrapper for ML models
library(UBL) # SMOTE

# visualisation
library(ggplot2) 
library(ggforce)
library(RColorBrewer)

##################################
rm(list=ls())

gc(reset = TRUE)
gc(reset = TRUE)
##################################
# Load a script
source("./Script/function/functions_util.R")
source("./Script/function/functions_plsr.R")

# ------------------------------------------------------------------------------
# Load study area shapefile and transform to UTM Zone 34S (EPSG:32734)
study_area <- st_read("./SHP/StudyArea2024.shp") %>% 
  st_transform(crs = 32734)

# ------------------------------------------------------------------------------
# Load ground-level photo metadata and convert to spatial data frame

# Note that "ground_level_photo_extracted.csv" is not included in the repository
# as it includes GPS information
sf_ground_photo <- read.csv("./CSV/ground_level_photo_extracted.csv") %>%
  mutate(
    # Parse JSON from .geo column to extract coordinates
    geo_parsed = map(.geo, ~ fromJSON(.)$coordinates),
    longitude = map_dbl(geo_parsed, 1),
    latitude = map_dbl(geo_parsed, 2)
  ) %>%
  dplyr::select(ID, longitude, latitude, BioRegion, Biome, Year, Month, 
                BareGround, all_of(vi),all_of(minerals),all_of(soil),all_of(topography)) %>%  # Drop unnecessary columns
  na.omit() %>%  # Remove rows with missing values (resulting in 990 rows)
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%  # Convert to sf object in WGS84
  st_transform(crs = 32734)  # Transform to UTM for spatial blocking

# ------------------------------------------------------------------------------
# Split data by growing season
sf_growing     <- sf_ground_photo %>% filter(Month < 5)
sf_non_growing <- sf_ground_photo %>% filter(Month >= 5)

# Selected Vegetation Indices to Apply Box-Cox transformation
growing_transform_vi <- c("bsi","evi","exg","ndmi","ndvi","swirRatio","tc_bright","tc_wet")
non_growing_transform_vi <- c("exg", "swirRatio","tc_wet")

# Get parameters for box-cox transformation for each VI during the growing season
growing_params <- map(growing_transform_vi, function(col) {
  params <- get_boxcox_params(sf_growing[[col]], sf_growing$BareGround)
  tibble(variable = col, lambda = params$lambda, shift_x = params$shift_x)
}) %>% 
  list_rbind()

# add columns for box-cox transformed values
sf_growing <- sf_growing %>%
  mutate(across(all_of(growing_transform_vi), 
                .fns = function(col) {
                  p <- growing_params %>% filter(variable == as.character(substitute(col)))
                  apply_boxcox_transform(col, p$lambda, p$shift_x)
                }, 
                .names = "{.col}_boxcox" 
                ))

# Get parameters for box-cox transformation for each VI during the non-growing season
non_growing_params <- map(non_growing_transform_vi, function(col) {
  params <- get_boxcox_params(sf_non_growing[[col]], sf_non_growing$BareGround)
  tibble(variable = col, lambda = params$lambda, shift_x = params$shift_x)
}) %>% 
  list_rbind()

# add columns for box-cox transformed values
sf_non_growing <- sf_non_growing %>%
  mutate(across(all_of(non_growing_transform_vi), 
                .fns = function(col) {
                  p <- non_growing_params %>% filter(variable == as.character(substitute(col)))
                  apply_boxcox_transform(col, p$lambda, p$shift_x)
                }, 
                .names = "{.col}_boxcox" 
  ))

# export box-cox parameter
write.csv(growing_params, "./Script/tmp/box-cox_growing_params.csv")
write.csv(non_growing_params, "./Script/tmp/box-cox_non_growing_params.csv")
rm(growing_params, non_growing_params)

# ------------------------------------------------------------------------------
# Data Split and export
sf_growing_split <- func_data_split(sf_growing, 0.8, k=5, seed=1234) #80:20 split & spatial blocking for CV
growing_tr_vi <- c(setdiff(vi, growing_transform_vi), paste0(growing_transform_vi, "_boxcox"))

sf_non_growing_split <- func_data_split(sf_non_growing, 0.8, k=5, seed=1234) #80:20 split & spatial blocking for CV
non_growing_tr_vi <- c(setdiff(vi, non_growing_transform_vi), paste0(non_growing_transform_vi, "_boxcox"))

saveRDS(sf_growing_split, file = "./Script/tmp/sf_growing_split.rds")
saveRDS(sf_non_growing_split, file = "./Script/tmp/sf_non_growing_split.rds")

# ------------------------------------------------------------------------------
# Run the pls model pipeline
library(plsVarSel)

#### 1-1. growing season #########################################################
#### PLS1-1 growing season without transformation #################################
full_growing <- caret::train(
  x=sf_growing_split$train_df[,vi], 
  y=sf_growing_split$train_df$BareGround,
  method="pls",
  preProcess = c('center', 'scale'),        # Standardize predictors before training
  importance = TRUE,                        # Enable variable importance calculation
  metric = "RMSE",
  trControl = trainControl(
    method = "cv",                          # Cross-validation
    selectionFunction = "oneSE",            # Select simplest model within 1-SE of best
    index = sf_growing_split$spatial_cv_block$indx_train, 
    indexOut = sf_growing_split$spatial_cv_block$indx_test,             # Cross-validation folds (testing)
    savePredictions = "final"               # Save final predictions
  ))

growing_results <- func_pls_var_selection(sf_growing_split, 
                                          vi, 
                                          full_growing$bestTune$ncomp)

#### PLS1-2 growing season with transformation ####################################
full_growing_tr <- caret::train(
  x=sf_growing_split$train_df[,growing_tr_vi], 
  y=sf_growing_split$train_df$BareGround,
  method="pls",
  preProcess = c('center', 'scale'),        # Standardize predictors before training
  importance = TRUE,                        # Enable variable importance calculation
  metric = "RMSE",
  trControl = trainControl(
    method = "cv",                          # Cross-validation
    selectionFunction = "oneSE",            # Select simplest model within 1-SE of best
    index = sf_growing_split$spatial_cv_block$indx_train, 
    indexOut = sf_growing_split$spatial_cv_block$indx_test,             # Cross-validation folds (testing)
    savePredictions = "final"               # Save final predictions
  ))

growing_tr_results <- func_pls_var_selection(sf_growing_split, 
                                          growing_tr_vi, 
                                          full_growing$bestTune$ncomp)

#### PLS1-3 growing season with SMOTE #############################################
full_growing_smote <- func_pls_smote(sf_growing_split$train_df, vi, sf_growing_split, seed=1234)
growing_smote_results <- func_pls_smote_var_selection(list_data=sf_growing_split,
                             predictors=vi,
                             n_comp=full_growing_smote$ncomp, 
                             seed=1234)

#### PLS1-4 growing season with SMOTE and transformation###########################
full_growing_tr_smote <- func_pls_smote(sf_growing_split$train_df, growing_tr_vi, sf_growing_split, seed=1234)
growing_tr_smote_results <- func_pls_smote_var_selection(list_data=sf_growing_split,
                                                      predictors=growing_tr_vi,
                                                      n_comp=full_growing_tr_smote$ncomp, 
                                                      seed=1234)

#### 1-2. non-growing season #####################################################
#### PLS2-1 non-growing season without transformation #############################
# run a full model
full_non_growing <- caret::train(
                          x=sf_non_growing_split$train_df[,vi], 
                          y=sf_non_growing_split$train_df$BareGround,
                          method="pls",
                          preProcess = c('center', 'scale'),        # Standardize predictors before training
                          importance = TRUE,                        # Enable variable importance calculation
                          metric = "RMSE",
                          trControl = trainControl(
                            method = "cv",                          # Cross-validation
                            selectionFunction = "oneSE",            # Select simplest model within 1-SE of best
                            index = sf_non_growing_split$spatial_cv_block$indx_train, 
                            indexOut = sf_non_growing_split$spatial_cv_block$indx_test,             # Cross-validation folds (testing)
                            savePredictions = "final"               # Save final predictions
                          ))
# full_non_growing
non_growing_results <- func_pls_var_selection(sf_non_growing_split, 
                                              vi, 
                                              full_non_growing$bestTune$ncomp)

#### PLS2-2 non-growing season with Box-Cox transformation ########################
# run a full model
full_non_growing_tr <- caret::train(
  x=sf_non_growing_split$train_df[,non_growing_tr_vi], 
  y=sf_non_growing_split$train_df$BareGround,
  method="pls",
  preProcess = c('center', 'scale'),        # Standardize predictors before training
  importance = TRUE,                        # Enable variable importance calculation
  metric = "RMSE",
  trControl = trainControl(
    method = "cv",                          # Cross-validation
    selectionFunction = "oneSE",            # Select simplest model within 1-SE of best
    index = sf_non_growing_split$spatial_cv_block$indx_train, 
    indexOut = sf_non_growing_split$spatial_cv_block$indx_test,             # Cross-validation folds (testing)
    savePredictions = "final"               # Save final predictions
  ))

non_growing_tr_results <- func_pls_var_selection(sf_non_growing_split, 
                                                 non_growing_tr_vi, 
                                                 full_non_growing$bestTune$ncomp)

# ------------------------------------------------------------------------------
# Results 

# ------------------------------------------------------------------------------
# Identify the best set of variables in different models
func_pls_select_var <- function(df_results){
  temp_df <- df_results$metrics %>% 
    slice(0:sum(df_results$vip_full_model<=1.0)+1) %>% # number of variables VIP >1.0
    mutate(n_var=(length(vi)-row_number()+1))

  min_idx   <- which.min(temp_df$RMSE)
  min_rmse  <- temp_df$RMSE[min_idx]
  se_rmse   <- temp_df$RMSESD[min_idx]/sqrt(5) # 5 fold cv
  selected_vars <- temp_df[temp_df$RMSE <= (min_rmse + se_rmse),] %>% 
    arrange(n_var)  %>% 
    dplyr::select(variables) %>% 
    head(1) %>% 
    unlist() %>% 
    unname()

  return(selected_vars)
}  

pls_growing_var          <- func_pls_select_var(growing_results)
pls_growing_tr_var       <- func_pls_select_var(growing_tr_results)
pls_growing_smote_var    <- func_pls_select_var(growing_smote_results)
pls_growing_tr_smote_var <- func_pls_select_var(growing_tr_smote_results)

pls_non_growing_var      <- func_pls_select_var(non_growing_results)
pls_non_growing_tr_var   <- func_pls_select_var(non_growing_tr_results)


# ------------------------------------------------------------------------------
# Run the best models again and save them
func_smote_final <- function(dataframe, predictors, blocks, n_comp){
  train_data <- dataframe[, c("BareGround",predictors)]
  smoter_train <- SmoteRegress(BareGround ~ ., dat = train_data, 
                               rel = create_relevance_matrix(train_data$BareGround), 
                               C.perc = "extreme",
                               k = 5)
  pls_model <- plsr(BareGround ~ ., data = smoter_train, ncomp=2, scale = TRUE, center = TRUE)
  return(pls_model)
}



pls_growing <- func_caret_pls(dataframe=sf_growing_split$train_df,
                               predictors = str_split(pls_growing_var, pattern = ",\\s*")[[1]],
                               blocks=sf_growing_split$spatial_cv_block,
                               n_comp=full_growing$bestTune$ncomp)

pls_growing_tr <- func_caret_pls(dataframe=sf_growing_split$train_df,
                              predictors = str_split(pls_growing_tr_var, pattern = ",\\s*")[[1]],
                              blocks=sf_growing_split$spatial_cv_block,
                              n_comp=full_growing_tr$bestTune$ncomp)

pls_growing_smote <- func_smote_final(dataframe=sf_growing_split$train_df,
                                 predictors = str_split(pls_growing_smote_var, pattern = ",\\s*")[[1]],
                                 blocks=sf_growing_split,
                                 n_comp=2)

pls_growing_tr_smote <- func_smote_final(dataframe=sf_growing_split$train_df,
                                      predictors = str_split(pls_growing_tr_smote_var, pattern = ",\\s*")[[1]],
                                      blocks=sf_growing_split,
                                      n_comp=2)

pls_non_growing <- func_caret_pls(dataframe=sf_non_growing_split$train_df,
                                   predictors = str_split(pls_non_growing_var, pattern = ",\\s*")[[1]],
                                   blocks=sf_non_growing_split$spatial_cv_block,
                                   n_comp=full_non_growing$bestTune$ncomp)

pls_non_growing_tr <- func_caret_pls(dataframe=sf_non_growing_split$train_df,
                                  predictors = str_split(pls_non_growing_tr_var, pattern = ",\\s*")[[1]],
                                  blocks=sf_non_growing_split$spatial_cv_block,
                                  n_comp=full_non_growing_tr$bestTune$ncomp)

# save models
saveRDS(pls_growing, file = "./Script/models/growing/pls_1-1_growing.rds")
saveRDS(pls_growing_tr, file = "./Script/models/growing/pls_1-2_growing_tr.rds")
saveRDS(pls_growing_smote, file = "./Script/models/growing/pls_1-3_growing_smote.rds")
saveRDS(pls_growing_tr_smote, file = "./Script/models/growing/pls_1-4_growing_tr_smote.rds")

saveRDS(pls_non_growing, file = "./Script/models/non_growing/pls_2-1_non_growing.rds")
saveRDS(pls_non_growing_tr, file = "./Script/models/non_growing/pls_2-2_non_growing_tr.rds")

# save predictors
list_predvars <- list(pls_growing_var,
                      pls_growing_tr_var,
                      pls_growing_smote_var,
                      pls_growing_tr_smote_var,
                      pls_non_growing_var,
                      pls_non_growing_tr_var)
out_list <- list()
for(predvars in list_predvars){
  temp <- paste0(predvars, collapse=", ")
  out_list <- c(out_list, temp)
}
summary_pls_model <- data.frame(
  model = paste0("pls_", c(rep(1,4),rep(2,2)), "-", c(1:4,1:2)),
  season = c(rep("growing",4),rep("non-growing",2)),
  preprocess = c("none", "box-cox", "smote", "box-cox & smote", "none", "box-cox"),
  ncomp=c(full_growing$bestTune$ncomp,
          full_growing_tr$bestTune$ncomp,
          full_growing_smote$ncomp,
          full_growing_tr_smote$ncomp,
          full_non_growing$bestTune$ncomp,
          full_non_growing_tr$bestTune$ncomp),
  predictors=unlist(out_list))

write.csv(summary_pls_model, "./Script/tmp/pls_summary.csv")

rm(full_growing,
   full_growing_tr,
   full_growing_smote,
   full_growing_tr_smote,
   full_non_growing,
   full_non_growing_tr,
   growing_results,
   growing_tr_results,
   growing_smote_results,
   growing_tr_smote_results,
   non_growing_results,
   non_growing_tr_results,
   list_predvars, 
   pls_growing_var,
   pls_growing_tr_var,
   pls_growing_smote_var,
   pls_growing_tr_smote_var,
   pls_non_growing_var,
   pls_non_growing_tr_var,
   out_list)
