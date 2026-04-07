# Script for 
#   - Table 4
#   - Table S3
#   - Figure 2
#   - Figure S2
#   - Table S4 and S5: ANCOVA to evaluate the topographic effects


# Load required libraries
library(tidyverse)
library(ggplot2)
library(car)

##################################
rm(list=ls())

gc(reset = TRUE)
gc(reset = TRUE)
##################################
# Load a script
source("./Script/function/functions_util.R")
source("./Script/function/functions_validation.R")

#-------------------------------------------------------------------------------
# Load models  
growing_models <- list.files("./Script/models/growing",
                             pattern = "\\.rds$", full.names=TRUE)

non_growing_models <- list.files("./Script/models/non_growing",
                                 pattern = "\\.rds$", full.names=TRUE)

#-------------------------------------------------------------------
# Photo validation dataset
df_photo_growing     <- readRDS("./Script/tmp/sf_growing_split.rds")$valid_df
df_photo_non_growing <- readRDS("./Script/tmp/sf_non_growing_split.rds")$valid_df

df_ph_gr_valid <- func_pred_valid(growing_models, df_photo_growing)
df_ph_ng_valid <- func_pred_valid(non_growing_models, df_photo_non_growing)

#----------------------------------------------------------------------
## Transect Data Validation
df_transect <- read.csv("./CSV/transect_extracted.csv")
growing_params <- read.csv("./Script/tmp/box-cox_growing_params.csv")
non_growing_params <- read.csv("./Script/tmp/box-cox_non_growing_params.csv")

# box-cox transformation for the selected variables
df_transect_growing <- df_transect %>%
  filter(Month < 5) %>%
  dplyr::select(BareGround, all_of(c(vi, soil,minerals,topography))) %>%
  mutate(across(all_of(growing_params$variable),
                .fns = function(col) {
                  p <- growing_params %>% filter(variable == as.character(substitute(col)))
                  apply_boxcox_transform(col, p$lambda, p$shift_x)
                },
                .names = "{.col}_boxcox"
  ))

df_transect_non_growing <- df_transect %>%
  filter(Month >= 5) %>%
  dplyr::select(BareGround, all_of(c(vi, soil,minerals,topography))) %>%
  mutate(across(all_of(non_growing_params$variable), 
                .fns = function(col) {
                  p <- non_growing_params %>% filter(variable == as.character(substitute(col)))
                  apply_boxcox_transform(col, p$lambda, p$shift_x)
                },
                .names = "{.col}_boxcox"
  ))

# Prediction using models 
df_tr_gr_valid <- func_pred_valid(growing_models, df_transect_growing)
df_tr_ng_valid <- func_pred_valid(non_growing_models, df_transect_non_growing)

#-------------------------------------------------------------------------------
# Table 4 Performance Metrics. R2, RMSE, MAE

# calculate evaluation metrics
photo_validation_summary <- rbind(
  func_eval_metrics(df_ph_gr_valid),
  func_eval_metrics(df_ph_ng_valid)
)
transect_validation_summary <- rbind(
  func_eval_metrics(df_tr_gr_valid),
  func_eval_metrics(df_tr_ng_valid)
)

# Format in a data frame
df_pls_summary      <- read.csv("./Script/tmp/pls_summary.csv")
df_pls_beta_summary <- read.csv("./Script/tmp/pls-beta_summary.csv")
df_rf_summary       <- read.csv("./Script/tmp/rf_summary.csv") %>% 
  mutate(X=row_number()) %>% 
  dplyr::select(X,model,season,preprocess,mtry,predictors)

photo_validation_summary <- bind_rows(
  photo_validation_summary %>% 
    filter(str_detect(model_name, "^pls_")) %>% 
    cbind(df_pls_summary),
  photo_validation_summary %>% 
    filter(str_detect(model_name, "^plsbeta_")) %>% 
    cbind(df_pls_beta_summary),
  photo_validation_summary %>% 
    filter(str_detect(model_name, "^rf_")) %>% 
    cbind(df_rf_summary)
  ) %>% 
  mutate(n_predictors = str_count(predictors, ",") + 1) 

transect_validation_summary <- bind_rows(
  transect_validation_summary %>% 
    filter(str_detect(model_name, "^pls_")) %>% 
    cbind(df_pls_summary),
  transect_validation_summary %>% 
    filter(str_detect(model_name, "^plsbeta_")) %>% 
    cbind(df_pls_beta_summary),
  transect_validation_summary %>% 
    filter(str_detect(model_name, "^rf_")) %>% 
    cbind(df_rf_summary)) %>% 
  mutate(n_predictors = str_count(predictors, ",") + 1) 

# write.csv(photo_validation_summary, "./Script/tmp/photo_validation.csv")
# write.csv(transect_validation_summary, "./Script/tmp/transect_validation.csv")

table4_photo <- photo_validation_summary %>% 
  mutate(Algorithm= str_split_i(model_name, "_",1)) %>% 
  dplyr::select(season, model_name,Algorithm,preprocess,n_predictors,rmse,mae,r_sq)

table4_transect <- transect_validation_summary %>% 
  mutate(Algorithm= str_split_i(model_name, "_",1)) %>% 
  dplyr::select(model_name,rmse,mae,r_sq)

table4 <- left_join(table4_photo, table4_transect, by="model_name") %>% 
  rename(photo_rmse=rmse.x,
         photo_mae=mae.x,
         photo_r_sq=r_sq.x,
         transect_rmse=rmse.y,
         transect_mae=mae.y,
         transect_r_sq=r_sq.y) %>% 
  arrange(season)

write.csv(table4, "./table/table4_performance_metrics.csv")

#-------------------------------------------------------------------------------
# Statistical test for evaluating the effect of preprocessing
# Table S3 

df_ph_lm <- table4 %>% 
  dplyr::select(-c(transect_rmse,transect_mae,transect_r_sq)) %>% 
  mutate(smote   = if_else(str_detect(preprocess, "smote"), 0, 1),
         box_cox = if_else(str_detect(preprocess, "box-cox"), 0, 1),
         season  = if_else(str_detect(model_name, "_1"), "non-growing", "growing")) %>% 
  mutate(Algorithm = as.factor(Algorithm),
         smote = as.factor(smote),
         box_cox = as.factor(box_cox),
         season = as.factor(season)) %>% 
  rename(rmse = photo_rmse,
         mae  = photo_mae,
         r_sq = photo_r_sq)

df_tr_lm <- table4 %>% 
  dplyr::select(-c(photo_rmse, photo_mae, photo_r_sq)) %>% 
  mutate(smote   = if_else(str_detect(preprocess, "smote"), 0, 1),
         box_cox = if_else(str_detect(preprocess, "box-cox"), 0, 1),
         season  = if_else(str_detect(model_name, "_1"), "non-growing", "growing")) %>% 
  mutate(Algorithm = as.factor(Algorithm),
         smote = as.factor(smote),
         box_cox = as.factor(box_cox),
         season = as.factor(season)) %>% 
  rename(rmse = transect_rmse,
         mae  = transect_mae,
         r_sq = transect_r_sq)

# anova test
ph_lm_table <- func_model_comparison(df_ph_lm) %>% 
  mutate(dataset = "Ground-level photograph")
tr_lm_table <- func_model_comparison(df_tr_lm) %>% 
  mutate(dataset = "Field-measured canopy cover")

# format the output tables
table_s3 <- rbind(ph_lm_table, tr_lm_table) 

# Export Table S3
write.csv(table_s3, "./table/table_s3_anova_output.csv")

#-------------------------------------------------------------------------------
# Draw scatter plots to evaluate model accuracy
# Figure 2 

  # -- Export settings --
png("./figures/figure2_scatter_plot_accuracy.png",
    width = 2400, height = 2200, res = 300)

par(
  mfrow = c(2, 2),
  mar   = c(2, 2, 2, 2),   # inner margins
  oma   = c(2, 4, 1, 0)    # outer margins
)

  # -- Growing season --

## Photo-based validation
df_ph_gr <- df_ph_gr_valid %>%
  na.omit()
plot_obs_pred(df_ph_gr$BareGround, df_ph_gr$`plsbeta_1-3`, "", y_lim = c(0, 100))

## Transect-based validation
df_tr_gr <- df_transect %>%
  filter(Month < 5) %>%                     # growing season
  dplyr::select(Notes, BareGround) %>%
  cbind(df_tr_gr_valid) %>%
  dplyr::select(-3) %>%                     # remove duplicate column
  na.omit()

plot_slope_flat(df_tr_gr, df_tr_gr$`plsbeta_1-3`, "", y_lim = c(0, 100))

legend(
  "bottomright",
  legend = c("Plain", "Slope"),
  title  = expression(bold("Landform")),
  col    = NA,
  pt.bg  = c(
    rgb(0, 0.8, 0.4, 0.6),
    rgb(1, 0, 0, 0.7)
  ),
  pch    = c(21, 24),
  pt.cex = 1.5,
  cex    = 1,
  bty    = "n"
)

  # -- Non-growing season --

## Photo-based validation
df_ph_ng <- df_ph_ng_valid %>%
  na.omit()
plot_obs_pred(df_ph_ng$BareGround, df_ph_ng$`plsbeta_2-1`, "", y_lim = c(0, 100))

## Transect-based validation
df_tr_ng <- df_transect %>%
  filter(Month > 4) %>%                     # non-growing season
  dplyr::select(Notes, BareGround) %>%
  cbind(df_tr_ng_valid) %>%
  dplyr::select(-3) %>%
  na.omit()

plot_slope_flat(df_tr_ng, df_tr_ng$`plsbeta_2-1`, "", y_lim = c(0, 100))

legend(
  "bottomright",
  legend = c("Plain", "Slope", "Others"),
  title  = expression(bold("Landform")),
  col    = NA,
  pt.bg  = c(
    rgb(0, 0.8, 0.4, 0.6),
    rgb(1, 0, 0, 0.7),
    rgb(0, 0, 0, 0.7)
  ),
  pch    = c(21, 24, 22),
  pt.cex = 1.5,
  cex    = 1,
  bty    = "n"
)

  # -- Common axis labels --
mtext("Predicted fraction of bare ground (%)", side  = 2, line = 0.5, outer = TRUE)
mtext("Observed fraction of bare ground (%)", side  = 1, line = 0.5, outer = TRUE)

  # -- Panel labels --
mtext("Ground-level Photograph", side = 3, line = -1, adj = 0.18, outer = TRUE)
mtext("Field-measured Canopy Cover", side = 3, line = -1, adj = 0.85, outer = TRUE)
mtext("Growing season", side = 2, line = 2.5, adj = 0.8, outer = TRUE)
mtext("Non-growing season", side = 2, line = 2.5, adj = 0.2, outer = TRUE)

  # -- Finish --

dev.off()

rm(df_ph_gr, df_tr_gr, df_ph_ng, df_tr_ng)
 
    
# ------------------------------------------------------------------------------
# Figure S2 

# Metadata for model configurations
# Algorithm names (PLS / PLS-beta / RF)
vec_algo <- rep(c("PLS", "PLS-beta", "RF"), 6, each = 2)

# Predictor sets used in each model
vec_predictors <- c(
  # non-growing
  rep("VI", 10), rep("VI & env vars", 2),
  rep("VI", 10), rep("VI & env vars", 2),
  # growing
  rep("VI", 10), rep("VI & env vars", 2)
)

# Pre-processing methods applied to predictors
vec_preprocess <- c(
  # non-growing
  rep("None", 6),
  rep("Box-Cox", 4), rep("None", 2),
  rep("SMOTE", 6),   
  rep("Box-Cox & SMOTE", 4), rep("SMOTE", 2),
  # growing
  rep("None", 6),
  rep("Box-Cox", 4), rep("None", 2)
)

# Model identifiers (used as panel titles)
vec_model_name <- c(
  paste0(vec_algo[1:24],  "1-", rep(1:4, each = 6)),
  paste0(vec_algo[25:36], "2-", rep(1:2, each = 6))
)

# Drawing order for panels
order_gr <- c(1, 5, 9,   # Order of models for growing season
              2, 6, 10,
              3, 7, 11,
              4, 8, 12)

order_ng <- c(1, 3, 5,  # Order of models for non-growing season
              2, 4, 6)

# ----- Setting for an export --------
png("./figures/figure_s2_scatter_plot_accuracy.png", width = 3200, height = 2800, res=200)  
# pdf(file = "./figures/figure_s2_scatter_plot_accuracy.pdf", width = 16, height = 14)

# ----- Drawing plots --------
# Layout settings for multi-panel plotting
par(
  mfrow = c(6, 6),          # 6 x 6 panel layout
  mar   = c(2, 2, 2, 2),    # inner margins
  oma   = c(2, 2, 0, 0)     # outer margins
)

# Growing season plots
for (i in 1:12) {
  
  cl <- order_gr[i] + 1     # column index for prediction
  j  <- i * 2 - 1           # index for model metadata
  
  ## ---- Photo-based validation --------------------------------
  df_ph_gr <- df_ph_gr_valid %>%
    na.omit()
  
  plot_obs_pred(
    df_ph_gr$BareGround,
    df_ph_gr[, cl],
    ""
  )
  
  mtext(paste0("algorithm: ",  vec_algo[j]),
        side = 1, adj = 0.95, line = -3.0, cex = 0.7)
  mtext(paste0("predictors: ", vec_predictors[j]),
        side = 1, adj = 0.95, line = -2.1, cex = 0.7)
  mtext(paste0("preprocess: ", vec_preprocess[j]),
        side = 1, adj = 0.95, line = -1.2, cex = 0.7)
  mtext(vec_model_name[j],
        side = 1, adj = 0.95, line = -4.3, cex = 0.9)
  
  ## ---- Transect-based validation ------------------------------
  df_tr_gr<- df_transect %>%
    filter(Month < 5) %>%                     # growing season
    dplyr::select(Notes, BareGround) %>%
    cbind(df_tr_gr_valid) %>%
    dplyr::select(-3) %>%                     # remove duplicate column
    na.omit()
  
  plot_slope_flat(
    df_tr_gr,
    df_tr_gr[, cl + 1],
    "transect-est"
  )
  
  mtext(paste0("algorithm: ",  vec_algo[j + 1]),
        side = 1, adj = 0.95, line = -3.0, cex = 0.7)
  mtext(paste0("predictors: ", vec_predictors[j + 1]),
        side = 1, adj = 0.95, line = -2.1, cex = 0.7)
  mtext(paste0("preprocess: ", vec_preprocess[j + 1]),
        side = 1, adj = 0.95, line = -1.2, cex = 0.7)
  mtext(vec_model_name[j + 1],
        side = 1, adj = 0.95, line = -4.3, cex = 0.9)
}


# Non-growing season plots
for (i in 1:6) {
  
  cl <- order_ng[i] + 1
  j  <- i * 2 - 1 + 24
  
  ## ---- Photo-based validation --------------------------------
  df_ph_ng <- df_ph_ng_valid %>%
    na.omit()
  
  plot_obs_pred(
    df_ph_ng$BareGround,
    df_ph_ng[, cl],
    ""
  )
  
  mtext(paste0("algorithm: ",  vec_algo[j]),
        side = 1, adj = 0.95, line = -3.0, cex = 0.7)
  mtext(paste0("predictors: ", vec_predictors[j]),
        side = 1, adj = 0.95, line = -2.1, cex = 0.7)
  mtext(paste0("preprocess: ", vec_preprocess[j]),
        side = 1, adj = 0.95, line = -1.2, cex = 0.7)
  mtext(vec_model_name[j],
        side = 1, adj = 0.95, line = -4.3, cex = 0.9)
  
  ## ---- Transect-based validation ------------------------------
  df_tr_ng <- df_transect %>%
    filter(Month > 4) %>%                     # non-growing season
    dplyr::select(Notes, BareGround) %>%
    cbind(df_tr_ng_valid) %>%
    dplyr::select(-3) %>%
    na.omit()
  
  plot_slope_flat(
    df_tr_ng,
    df_tr_ng[, cl + 1],
    "transect-est"
  )
  
  mtext(paste0("algorithm: ",  vec_algo[j + 1]),
        side = 1, adj = 0.95, line = -3.0, cex = 0.7)
  mtext(paste0("predictors: ", vec_predictors[j + 1]),
        side = 1, adj = 0.95, line = -2.1, cex = 0.7)
  mtext(paste0("preprocess: ", vec_preprocess[j + 1]),
        side = 1, adj = 0.95, line = -1.2, cex = 0.7)
  mtext(vec_model_name[j + 1],
        side = 1, adj = 0.95, line = -4.3, cex = 0.9)
}


# Common axis labels (outer margins)
mtext(
  "Predicted fraction of bare ground (%)",
  side  = 2,
  line  = 0.5,
  outer = TRUE
)

mtext(
  "Observed fraction of bare ground (%)",
  side  = 1,
  line  = 0.5,
  outer = TRUE
)
dev.off() 

rm(df_ph_gr, df_tr_gr, df_ph_ng, df_tr_ng)


# ------------------------------------------------------------------------------
# Pariwise Tukey test, 
# Figure S3

df_gg_posthoc_ph <- func_posthoc_tukey(df_ph_lm) %>% 
  mutate(label = paste(pairs, metrics, sep=": "),
         x_pos = c(1.5,1.25,1, 3,2.75,2.5, 4.5,4.25,4)) 

df_gg_posthoc_tr <- func_posthoc_tukey(df_tr_lm) %>% 
  mutate(label = paste(pairs, metrics, sep=": "),
         x_pos = c(1.5,1.25,1, 3,2.75,2.5, 4.5,4.25,4)) 

ph_tukey_test <- gg_posthoc(df_gg_posthoc_ph)
tr_tukey_test <- gg_posthoc(df_gg_posthoc_tr)

# Export plot
ggsave(file="./figures/figure_s3a_tukey_photo.png",
       plot=ph_tukey_test,
       device="png",
       width=2500, height=1200,
       units="px")

ggsave(file="./figures/figure_s3b_tukey_transect.png",
       plot=tr_tukey_test,
       device="png",
       width=2500, height=1200,
       units="px")


#-------------------------------------------------------------------------------
# Table S4 Statistical Test to see the effects of slope vs flat plains
# ANCOVA

colnames(df_tr_gr_valid) <- gsub("-", "_", colnames(df_tr_gr_valid))
colnames(df_tr_ng_valid) <- gsub("-", "_", colnames(df_tr_ng_valid))

df_gr_ancova <- df_transect %>%
  filter(Month < 5) %>% 
  dplyr::select(Notes, BareGround) %>% 
  cbind(df_tr_gr_valid) %>% 
  dplyr::select(-3) %>% 
  mutate(landforms = as.factor(Notes)) %>% 
  filter(landforms !="Others")

df_ng_ancova <- df_transect %>%
  filter(Month > 4) %>% 
  dplyr::select(Notes, BareGround) %>% 
  cbind(df_tr_ng_valid) %>% 
  dplyr::select(-3) %>% 
  mutate(landforms = as.factor(Notes)) %>% 
  filter(landforms !="Others")

gr_ancova<- func_ancova(df_gr_ancova) 
ng_ancova<- func_ancova(df_ng_ancova) 

# export
write.csv(gr_ancova,"./table/table_s4_anova_topography_growing.csv")
write.csv(ng_ancova,"./table/table_s5_anova_topography_non-growing.csv")
