# This is a script to draw graphs for Figure S3

# Load required libraries
library(tidyverse)
library(ggplot2)
library(patchwork)
##################################
rm(list=ls())

gc(reset = TRUE)
gc(reset = TRUE)
##################################

# Load a script
source("./Script/function/functions_util.R")
source("./Script/function/functions_residual_analysis.R")

# Define predictor variable groups
env_cov <- c("Biome", "BioRegion","Year", "Month")
#-------------------------------------------------------------------------------
# Load models and datasets
plsbeta_1_3 <- readRDS("./Script/models/growing/plsbeta_1-3_growing_smote.rds")
plsbeta_2_1 <- readRDS("./Script/models/non_growing/plsbeta_2-1_non_growing.rds")

# Photo validation dataset
sf_growing_split <- readRDS("./Script/tmp/sf_growing_split.rds")
sf_non_growing_split <- readRDS("./Script/tmp/sf_non_growing_split.rds")

## Transect Data Validation
df_transect <- read.csv("./CSV/transect_extracted.csv")

#-------------------------------------------------------------------------------
# Data wrangling for ggplot

# Photo Dataset
# make a prediction using loaded data and models for photograph data
df_ph_gr_valid <- func_pred_valid(plsbeta_1_3, sf_growing_split$valid_df) %>% 
  mutate(data_type = "photo_valid")
df_ph_ng_valid <- func_pred_valid(plsbeta_2_1, sf_non_growing_split$valid_df) %>% 
  mutate(data_type = "photo_valid")

df_ph_gr_train <- func_pred_valid(plsbeta_1_3, sf_growing_split$train_df) %>% 
  mutate(data_type = "photo_train")
df_ph_ng_train <- func_pred_valid(plsbeta_2_1, sf_non_growing_split$train_df) %>% 
  mutate(data_type = "photo_train")

# Transect Dataset 
df_transect_growing <- df_transect %>%
  filter(Month < 5) %>%
  dplyr::select(BareGround, all_of(vi), all_of(env_cov)) 

df_transect_non_growing <- df_transect %>%
  filter(Month >= 5) %>%
  dplyr::select(BareGround, all_of(vi), all_of(env_cov)) 

df_tr_gr_valid <- func_pred_valid(plsbeta_1_3, df_transect_growing) %>% 
  mutate(data_type = "transect")
df_tr_ng_valid <- func_pred_valid(plsbeta_2_1, df_transect_non_growing) %>% 
  mutate(data_type = "transect")

# combine photo and transect dataset to make a dataframe for ggplot
df_gr_res <- rbind(df_ph_gr_train, df_ph_gr_valid, df_tr_gr_valid) %>% 
  mutate(residuals = BareGround-predicted)
df_ng_res <- rbind(df_ph_ng_train, df_ph_ng_valid, df_tr_ng_valid) %>% 
  mutate(residuals = BareGround-predicted)

#-------------------------------------------------------------------------------
# Drawing Figure S8 Residuals vs Covariate

# By year
gg_res_gr_year <- df_gr_res %>% 
  filter(data_type == "photo_train" | data_type == "photo_valid") %>% 
  plot_residuals_by_covariate(color_column = "Year", color_option = "D", color_label = "Year", color_lim=c(2014,2024)) +
  theme(legend.position = "none") +
  annotate("text", x=5, y=50, label="A", fontface="bold", size= 7) 

gg_res_ng_year <- df_ng_res %>% 
  filter(data_type == "photo_train" | data_type == "photo_valid") %>% 
  plot_residuals_by_covariate(color_column = "Year", color_option = "D", color_label = "Year", color_lim=c(2014,2024)) +
  theme(legend.position = "right") +
  annotate("text", x=5, y=50, label="B", fontface="bold", size= 7) 

# By month
gg_res_gr_month <- df_gr_res %>% 
  filter(data_type == "photo_train" | data_type == "photo_valid") %>% 
  plot_residuals_by_covariate(color_column = "Month", color_option = "C", color_label = "Month", color_lim=c(1,12))+
  theme(legend.position = "none") +
  annotate("text", x=5,   y=50, label="C", fontface="bold", size= 7) 

gg_res_ng_month <- df_ng_res %>% 
  filter(data_type == "photo_train" | data_type == "photo_valid") %>% 
  plot_residuals_by_covariate(color_column = "Month", color_option = "C", color_label = "Month", color_lim=c(1,12))+
  theme(legend.position = "right")+
  annotate("text", x=5, y=50, label="D", fontface="bold", size= 7) 


# By Biome and Bioregion
color_mapping <- c(
  "Bushmanland Bioregion" = "orange",
  "Lower Karoo Bioregion" = "brown",
  "Upper Karoo Bioregion" = "hotpink",
  "Desert" = "grey60",
  "Grassland" = "green2",
  "Albany Thicket" = "darkolivegreen",
  "Azonal Vegetation" = "blue1"
)

label_mapping <- c(
  "orange" = "Bushmanland",
  "brown" = "Lower Karoo",
  "hotpink" = "Upper Karoo",
  "grey60" = "Desert",
  "green2" = "Grassland",
  "darkolivegreen" = "Albany Thicket",
  "blue1" = "Azonal Vegetation"
)

gg_res_gr_biome <- df_gr_res %>%
  plot_residuals_by_biome()+
  theme(legend.position = "none")+
  annotate("text", x=5, y=50, label="E", fontface="bold", size= 7) 

gg_res_ng_biome <- df_ng_res %>%
  plot_residuals_by_biome()+
  theme(legend.position = "right",
        legend.key.size = unit(0.9, "lines")
        )+
  annotate("text", x=5, y=50, label="F", fontface="bold", size= 7) 

res_plot <- (gg_res_gr_year + gg_res_ng_year)/
  (gg_res_gr_month + gg_res_ng_month)/
  (gg_res_gr_biome +   gg_res_ng_biome)

# export graphs
ggsave(file="./figures/figure_s8_Residual_analysis.png",
       plot= res_plot,
       device="png",
       width=2800, height=3200,
       units="px")



rm(gg_res_gr_month, gg_res_ng_month,
   gg_res_gr_year,  gg_res_ng_year,
   gg_res_gr_biome, gg_res_ng_biome)

