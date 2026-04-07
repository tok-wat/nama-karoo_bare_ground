library(tidyverse)
library(UBL)

##################################
rm(list=ls())

gc(reset = TRUE)
gc(reset = TRUE)

##################################

# Load a script
source("./Script/function/functions_util.R")
source("./Script/function/functions_pls_beta.R")

# ------------------------------------------------------------------------------
# Load preprocessed file

sf_growing_split <- readRDS("./Script/tmp/sf_growing_split.rds")
sf_non_growing_split <- readRDS("./Script/tmp/sf_non_growing_split.rds")

# nrow(sf_growing_split$train_df)
# nrow(sf_non_growing_split$train_df)


vi <- c("bsi","evi","exg","msavi","ndmi","ndvi","osavi","satvi","savi","swirRatio","tc_bright","tc_green","tc_wet")
# Selected Vegetation Indices to Apply Box-Cox transformation
growing_transform_vi <- c("bsi","evi","exg","ndmi","ndvi","swirRatio","tc_bright","tc_wet")
non_growing_transform_vi <- c("exg", "swirRatio","tc_wet")

growing_tr_vi <- c(setdiff(vi, growing_transform_vi), paste0(growing_transform_vi, "_boxcox"))
non_growing_tr_vi <- c(setdiff(vi, non_growing_transform_vi), paste0(non_growing_transform_vi, "_boxcox"))


# ------------------------------------------------------------------------------
# Run the pls-beta model pipeline
library(plsRbeta)

# Preprocess training data for beta regression
# Scale BareGround to the range [0,1]
# Adjust 0 and 1 values with 1e-8 correction
df_plsbeta_growing <- sf_growing_split$train_df  %>% 
  mutate(BareGround = BareGround/100) %>% 
  mutate(BareGround = if_else(BareGround == 0, 0.001, BareGround)) %>% 
  mutate(BareGround = if_else(BareGround == 1, 0.999, BareGround))

df_plsbeta_non_growing <- sf_non_growing_split$train_df  %>% 
  mutate(BareGround = BareGround/100) %>% 
  mutate(BareGround = if_else(BareGround == 0, 0.001, BareGround)) %>% 
  mutate(BareGround = if_else(BareGround == 1, 0.999, BareGround))

#### 2-1. growing season #########################################################
#### PLS-beta1-1 growing season without transformation #################################

# Determine the optimal number of components (Nb_comp) from cross-validation results
plsbeta_growing_cv <- func_plsRbeta_ncomp_cv(
  dataframe = df_plsbeta_growing,
  predictors = vi,
  data_split = sf_growing_split)

# according to the help file, Small AIC & BIC, Large Q2Chisq_Y, Chi2_Pearson_Y, pseudo_R2_Y, R2_Y
# indicate good model but there is no consensus which metrics works the best
# https://fbertran.github.io/plsRbeta/articles/tn_Insights.html
# Here I used 'kink' method to select number of components noting AIC and BIC often
# overestimate numbers of components. Note that AIC, BIC, PseudoR2, and R2 use full
# model for calculation (not using CV). 

# par(mfrow=c(2,3))
# plot_cv(plsbeta_growing_cv)

# Perform bootstrapping for PLS-beta model with the selected number of components
growing_ncomp <- 3

plsbeta_growing_init <- plsRbeta(
  BareGround ~ ., 
  data = df_plsbeta_growing[,c("BareGround", vi)], 
  nt = growing_ncomp, 
  K=5, grouplist=list(sf_growing_split$spatial_cv_block$indx_test),
  modele = "pls-beta", verbose = FALSE) %>% 
  bootplsbeta(typeboot = "fmodel_np", R = 1000)

# Visualize bootstrap output for variable selection
# Select significant variables based on confidence intervals
# plsRglm::boxplots.bootpls(plsbeta_growing_init)
# plot_plsbeta_boot(plsbeta_growing_init)

plsbeta_growing_vars <- plsRglm::confints.bootpls(plsbeta_growing_init, indices = 1:length(vi)) %>% 
  get_selected_var() # Select using BCa confidence intervals

# Rebuild PLS-beta model using selected variables
plsbeta_growing <- plsRbeta(
  BareGround ~ .,
  data = df_plsbeta_growing[,c("BareGround", plsbeta_growing_vars)], 
  nt = growing_ncomp, 
  modele = "pls-beta",
  link.phi= "log",
  scaleX = TRUE,
  verbose = FALSE
)

#### PLS-beta1-2 growing season with transformation ####################################
# Determine the optimal number of components (Nb_comp) from cross-validation results
plsbeta_growing_tr_cv <- func_plsRbeta_ncomp_cv(
  dataframe = df_plsbeta_growing,
  predictors = growing_tr_vi,
  data_split = sf_growing_split)

# plot_cv(plsbeta_growing_tr_cv)

# Run Pls-beta regression
growing_tr_ncomp <- 3
plsbeta_growing_tr_init <- plsRbeta(
  BareGround ~ ., 
  data = df_plsbeta_growing[,c("BareGround", growing_tr_vi)], 
  nt = growing_tr_ncomp, 
  K=5, grouplist=list(sf_growing_split$spatial_cv_block$indx_test),
  modele = "pls-beta", verbose = FALSE) %>%
  bootplsbeta(typeboot = "fmodel_np", R = 1000)

# variable selection based on bootstrap
plsbeta_growing_tr_vars  <- plsRglm::confints.bootpls(plsbeta_growing_tr_init, indices = 1:length(growing_tr_vi)) %>% 
  get_selected_var() # Select using BCa confidence intervals

# Rebuild PLS-beta model using selected variables
plsbeta_growing_tr <- plsRbeta(
  BareGround ~ .,
  data = df_plsbeta_growing[,c("BareGround", plsbeta_growing_tr_vars)], 
  nt = growing_tr_ncomp, 
  modele = "pls-beta",
  link.phi= "log",
  scaleX = TRUE,
  K = 5, 
  grouplist= list(sf_growing_split$spatial_cv_block$indx_test),
  verbose = FALSE
)

#### PLS-beta1-3 growing season with SMOTE #############################################
# Determine the optimal number of components (Nb_comp) from cross-validation results
plsbeta_growing_smote_cv <- func_plsRbeta_smote_ncomp_cv(
  dataframe=df_plsbeta_growing,
  predictors=vi,
  data_split=sf_growing_split,
  seed=1234)

# plot_cv(plsbeta_growing_smote_cv)

# Run Pls-beta regression
growing_smote_ncomp <- 3
set.seed(1234)
df_growing_smoter <- SmoteRegress(
  BareGround ~ ., dat = df_plsbeta_growing[,c("BareGround",vi)], 
  rel = create_relevance_matrix(df_plsbeta_growing$BareGround), 
  C.perc = "extreme",
  k = 5)

plsbeta_growing_smote_init <- plsRbeta(
  BareGround ~ ., 
  data = df_growing_smoter[,c("BareGround", vi)], 
  nt = growing_smote_ncomp, 
  modele = "pls-beta", verbose = FALSE) %>%
  bootplsbeta(typeboot = "fmodel_np", R = 1000)

plsbeta_growing_smote_vars <- plsRglm::confints.bootpls(
  plsbeta_growing_smote_init, 
  indices = 1:length(vi)) %>% 
  get_selected_var() 

# Rebuild PLS-beta model using selected variables
plsbeta_growing_smote <- plsRbeta(
  BareGround ~ .,
  data = df_growing_smoter[,c("BareGround", plsbeta_growing_smote_vars)], 
  nt = growing_smote_ncomp, 
  modele = "pls-beta",
  link.phi= "log",
  scaleX = TRUE,
  verbose = FALSE
)

#### PLS-beta1-4 growing season with Box-Cox transformation and SMOTE ###################
# Determine the optimal number of components (Nb_comp) from cross-validation results
plsbeta_growing_tr_smote_cv <- func_plsRbeta_smote_ncomp_cv(
  dataframe=df_plsbeta_growing,
  predictors=growing_tr_vi,
  data_split=sf_growing_split,
  seed=1234)

# plot_cv(plsbeta_growing_tr_smote_cv)

# Run Pls-beta regression
growing_tr_smote_ncomp <- 3
set.seed(1234)
df_growing_tr_smote <- SmoteRegress(
  BareGround ~ ., dat = df_plsbeta_growing[,c("BareGround",growing_tr_vi)], 
  rel = create_relevance_matrix(df_plsbeta_growing$BareGround), 
  C.perc = "extreme",
  k = 5)

plsbeta_growing_tr_smote_init <- plsRbeta(
  BareGround ~ ., 
  data = df_growing_tr_smote[,c("BareGround", growing_tr_vi)], 
  nt = growing_tr_smote_ncomp, 
  modele = "pls-beta", verbose = FALSE) %>%
  bootplsbeta(typeboot = "fmodel_np", R = 1000)

plsbeta_growing_smote_tr_vars <- plsRglm::confints.bootpls(
  plsbeta_growing_tr_smote_init, 
  indices = 1:length(growing_tr_vi)) %>% 
  get_selected_var() 

# Rebuild PLS-beta model using selected variables
plsbeta_growing_tr_smote <- plsRbeta(
  BareGround ~ .,
  data = df_growing_tr_smote[,c("BareGround", plsbeta_growing_smote_tr_vars)], 
  nt = growing_tr_smote_ncomp, 
  modele = "pls-beta",
  link.phi= "log",
  scaleX = TRUE,
  verbose = FALSE
)

#### 2-2. non-growing season #####################################################
#### PLS-beta2-1 non-growing season without transformation #############################
plsbeta_non_growing_cv <- func_plsRbeta_ncomp_cv(
  dataframe = df_plsbeta_non_growing,
  predictors = vi,
  data_split = sf_non_growing_split,
  seed=1234)

# plot_cv(plsbeta_non_growing_cv)

# Run Pls-beta regression
non_growing_ncomp <- 2
plsbeta_non_growing_init <- plsRbeta(
  BareGround ~ ., data = df_plsbeta_non_growing[,c("BareGround", vi)], 
  nt = non_growing_ncomp, 
  K=5, grouplist=list(sf_non_growing_split$spatial_cv_block$indx_test),
  modele = "pls-beta", verbose = FALSE
) %>% bootplsbeta(typeboot = "fmodel_np", R = 1000)

# Select significant variables based on confidence intervals
plsbeta_non_growing_vars <- plsRglm::confints.bootpls(plsbeta_non_growing_init, indices = 1:length(vi)) %>% 
  get_selected_var() # Select using BCa confidence intervals

# Rebuild PLS-beta model using selected variables
plsbeta_non_growing <- plsRbeta(
  BareGround ~ .,
  data = df_plsbeta_non_growing[,c("BareGround", plsbeta_non_growing_vars)], 
  nt = non_growing_ncomp, 
  modele = "pls-beta",
  link.phi= "log",
  scaleX = TRUE,
  K = 5, 
  grouplist= list(sf_non_growing_split$spatial_cv_block$indx_test),
  verbose = FALSE
)

#### PLS-beta2-2 non-growing season with Box-Cox transformation ########################
plsbeta_non_growing_tr_cv <- func_plsRbeta_ncomp_cv(
  dataframe = df_plsbeta_non_growing,
  predictors = non_growing_tr_vi,
  data_split = sf_non_growing_split,
  seed=1234)

# plot_cv(plsbeta_non_growing_tr_cv)

# Run Pls-beta regression
non_growing_tr_ncomp <- 2
plsbeta_non_growing_tr_init <- plsRbeta(
  BareGround ~ ., data = df_plsbeta_non_growing[,c("BareGround", non_growing_tr_vi)], 
  nt = non_growing_tr_ncomp, 
  K=5, grouplist=list(sf_non_growing_split$spatial_cv_block$indx_test),
  modele = "pls-beta", verbose = FALSE
) %>% bootplsbeta(typeboot = "fmodel_np", R = 1000)


# Select significant variables based on confidence intervals
plsbeta_non_growing_tr_vars <- plsRglm::confints.bootpls(
  plsbeta_non_growing_tr_init, 
  indices = 1:length(non_growing_tr_vi)) %>% 
  get_selected_var() # Select using BCa confidence intervals

# Rebuild PLS-beta model using selected variables
plsbeta_non_growing_tr <- plsRbeta(
  BareGround ~ .,
  data = df_plsbeta_non_growing[,c("BareGround", plsbeta_non_growing_tr_vars)], 
  nt = non_growing_tr_ncomp, 
  modele = "pls-beta",
  link.phi= "log",
  scaleX = TRUE,
  K = 5, 
  grouplist= list(sf_non_growing_split$spatial_cv_block$indx_test),
  verbose = FALSE
)

# ------------------------------------------------------------------------------
# Export 

# save models
saveRDS(plsbeta_growing, file = "./Script/models/growing/plsbeta_1-1_growing.rds")
saveRDS(plsbeta_growing_tr, file = "./Script/models/growing/plsbeta_1-2_growing_tr.rds")
saveRDS(plsbeta_growing_smote, file = "./Script/models/growing/plsbeta_1-3_growing_smote.rds")
saveRDS(plsbeta_growing_tr_smote, file = "./Script/models/growing/plsbeta_1-4_growing_tr_smote.rds")

saveRDS(plsbeta_non_growing, file = "./Script/models/non_growing/plsbeta_2-1_non_growing.rds")
saveRDS(plsbeta_non_growing_tr, file = "./Script/models/non_growing/plsbeta_2-2_non_growing_tr.rds")

# save predictors
list_predvars <- list(plsbeta_growing_vars,
     plsbeta_growing_tr_vars,
     plsbeta_growing_smote_vars,
     plsbeta_growing_smote_tr_vars,
     plsbeta_non_growing_vars,
     plsbeta_non_growing_tr_vars)
out_list <- list()
for(predvars in list_predvars){
  temp <- paste0(predvars, collapse=", ")
  out_list <- c(out_list, temp)
}
summary_pls_model <- data.frame(
  model = paste0("plsbeta_", c(rep(1,4),rep(2,2)), "-", c(1:4,1:2)),
  season = c(rep("growing",4),rep("non-growing",2)),
  preprocess = c("none", "box-cox", "smote", "box-cox & smote", "none", "box-cox"),
  ncomp=c(growing_ncomp,growing_tr_ncomp,growing_smote_ncomp, growing_tr_smote_ncomp,
          non_growing_ncomp, non_growing_tr_ncomp),
  predictors=unlist(out_list))
write.csv(summary_pls_model, "./Script/tmp/pls-beta_summary.csv")

