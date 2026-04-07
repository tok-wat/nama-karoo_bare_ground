# Script for 
#   - Figure 3 Biplot
#   - Figure 4 Comparison with Google Earth Image
#   - Figure S1 Comparison of transect and photographic mesurement

library(tidyverse)
library(rlang)
library(jsonlite)

# spatial mapping
library(sf)
library(raster)
library(tidyterra)

# modelling
library(plsVarSel)
library(randomForest)

# visualisation
library(ggplot2)
library(grid)
library(RColorBrewer)
library(patchwork)
library(magick)
library(cowplot)
##################################
rm(list=ls())

gc(reset = TRUE)
gc(reset = TRUE)
##################################
# Load a script
source("./Script/function/functions_visualisation.R")

# ------------------------------------------------------------------------------
#  Figure 3 PLSR-beta Biplot
best_growing <- readRDS("./script/models/growing/plsbeta_1-3_growing_smote.rds")
best_non_growing <- readRDS("./script/models/non_growing/plsbeta_2-1_non_growing.rds")

gg_biplot_gr_12 <- plot_biplot(best_growing, label_scale=1.2, score_scale=1, components=1:2, show_legend=FALSE) +
  annotate("text", x=-5,   y=0.9, label="A",
           fontface="bold", size=5) 
gg_biplot_gr_13 <- plot_biplot(best_growing, 1.2, 1, c(1,3), show_legend=FALSE) +
  annotate("text", x=-5,   y=0.9, label="B",
           fontface="bold", size=5) 
gg_biplot_ng_12 <- plot_biplot(best_non_growing, 1.2, 7, 1:2) +
  annotate("text", x=-7,   y=6.3, label="D",
           fontface="bold", size=5) 

gg_combined_biplot <- gg_biplot_gr_12 + gg_biplot_ng_12 + gg_biplot_gr_13 +  plot_layout(ncol = 2,  guides = 'collect') 

# prepare standardised coefficient
df_std_coef_gr <- process_std_coef(best_growing$Std.Coeffs, vi_labels, vi_colours)
df_std_coef_ng <- process_std_coef(best_non_growing$Std.Coeffs, vi_labels, vi_colours)

# plot standardised coefficient
gg_std_coef_gr <- plot_std_coef(df_std_coef_gr) + 
  annotate("text", x=-0.8,  y= 4.3, label="C",  fontface="bold", size=5) 
gg_std_coef_ng <- plot_std_coef(df_std_coef_ng) +
  annotate("text", x=-0.8,  y= (nrow(df_std_coef_ng)+1)*0.9, label="E",  fontface="bold", size=5) 
gg_combined_coefplot <- gg_std_coef_gr | gg_std_coef_ng

# export plots
gg_combined_biplot_coefplot <- gg_combined_biplot/gg_combined_coefplot + plot_layout(heights=c(3,1))
ggsave(file="./figures/figure3_biplot_coefplot.pdf",
       plot=gg_combined_biplot_coefplot,
       device="pdf",
       width=2570, height=2600,
       units="px")

# ------------------------------------------------------------------------------
# Figure 4 Spatial Application

# Function for drawing three maps (rows) at a certain place
three_plot <- function(place){
  # load the shp file and then transform
  shp <- st_read("./SHP/Figure_4_Map_5km.shp")
  shp_utm <- st_transform(shp, crs = 32734)
  shp_utm <- shp_utm %>% filter(location==place)

  rgb_path <- paste0("./TIFF/google_earth_",place,".tif")
  val1_path <- paste0("./TIFF/pred_Mar_2024_",place,".tif")
  val2_path <- paste0("./TIFF/pred_Oct_2024_",place,".tif")
  
  df_rgb <- prepare_rgb_raster(rgb_path, shp_utm)
  df_val1 <- read_and_crop_raster(val1_path, shp_utm)
  df_val2 <- read_and_crop_raster(val2_path, shp_utm)
  
  p1 <- plot_ge(df_rgb, " ")
  p2 <- plot_rast(df_val1, " ")
  p3 <- plot_rast(df_val2, " ")
  
  return(
    list(p1=p1, p2=p2, p3=p3)
  )
}

# labelling and ggplot specifications
# Column labels
col1 <- ggplot() + ggtitle("VHR satellite images\n(Google Eatrh)") + theme_void() + 
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))
col2 <- ggplot() + ggtitle("March 2024\n(Growing season)") + theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))
col3 <- ggplot() + ggtitle("October 2024\n(Non-growing season)") + theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))
header_row <- wrap_elements(grid::nullGrob()) + col1 + col2 + col3 + 
  plot_layout(ncol = 4, guides = "collect", widths = c(0.1, 1, 1, 1)) & theme(plot.margin = ggplot2::margin(3, 3, 0.1, 3)) 

# Row labels
label_a <- wrap_elements(grid::textGrob("NW Bushmanland", rot = 90, gp = gpar(fontsize = 12, fontface = "bold")))
label_b <- wrap_elements(grid::textGrob("Carnarvon", rot = 90, gp = gpar(fontsize = 12, fontface = "bold")))
label_c <- wrap_elements(grid::textGrob("Leeu Gamka", rot = 90, gp = gpar(fontsize = 12, fontface = "bold")))

# making rows, adding row labels
row1 <- label_a + 
  three_plot("R358")$p1 + three_plot("R358")$p2 + three_plot("R358")$p3 + 
  plot_layout(ncol = 4, guides = "collect", widths = c(0.1, 1, 1, 1)) & theme(plot.margin = ggplot2::margin(3, 3, 3, 3)) 
row2 <- label_b + 
  three_plot("Carnarvon")$p1 + three_plot("Carnarvon")$p2 + three_plot("Carnarvon")$p3 +
  plot_layout(ncol = 4, guides = "collect", widths = c(0.1, 1, 1, 1)) & theme(plot.margin = ggplot2::margin(3, 3, 3, 3)) 
row3 <- label_c +  
  three_plot("LeeuGamka")$p1 + three_plot("LeeuGamka")$p2 + three_plot("LeeuGamka")$p3 + 
  plot_layout(ncol = 4, guides = "collect", widths = c(0.1, 1, 1, 1)) & theme(plot.margin = ggplot2::margin(3, 3, 3, 3)) 

#  make legend only
temp_gg <- data.frame(
  x = 1:10,
  y = 1:10,
  value = seq(0, 100, length.out = 10)) %>% 
  ggplot(aes(x = x, y = y, fill = value)) +
  geom_raster() +
  scale_fill_gradientn(
    colors = rev(brewer.pal(11, "BrBG")),
    limits = c(0, 100),
    oob = scales::squish,
    name = "Bare Ground (%)"
  ) +
  theme_void() +
  theme(legend.position = "right")  

legend <- cowplot::get_legend(temp_gg)
ggsave("./figures/figure4_legend_only.png", plot = legend, width = 2, height = 4, dpi = 300)

# Patchwork them
gg_comparison <- (header_row) / row1 / row2 / row3 +
  plot_layout(ncol = 1, heights = c(0.07, 1, 1, 1), widths = c(0.15, 1, 1, 1), guides="collect")

# export as pdf
ggsave(file="./figures/figure4_enlarged_prediction_map.pdf",
       plot=gg_comparison,
       device="pdf",
       width=1900, height=2000,
       units="px")

# ------------------------------------------------------------------------------
# Figure S1 Photo-transect comparison
df_transect_photo_comparison <- read.csv("./CSV/transect-photo_comparison.csv")

png("./figures/figure_s1_transect_photo_compare.png", width = 1000, height = 1100, res=200)  
  plot(df_transect_photo_comparison$ph_bare,
       df_transect_photo_comparison$tr_bare,
       xlab="Field measured fraction of bare ground (%)",
       ylab="Photo estimated fraction of bare ground (%)",
       xlim =c(0,100),ylim=c(0,100),
       col=NULL, bg=rgb(0, 0, 0, alpha=0.3), pch=21, cex=1.5
        )
  # decorate the plot
  temp_lm <-lm(df_transect_photo_comparison$tr_bare ~ df_transect_photo_comparison$ph_bare) # Rsq=0.79, p<0.001
  abline(0,1, col="black", lwd=2, lty=2)
  abline(temp_lm, col="red", lwd=2)
  r_sq <- round(summary(temp_lm)$r.squared,3)
  rmse <- mean((df_transect_photo_comparison$ph_bare-df_transect_photo_comparison$tr_bare)^2) %>%
    sqrt() %>%
    round(3)
  mae <- mean(abs(df_transect_photo_comparison$ph_bare-df_transect_photo_comparison$tr_bare)) %>%
    round(3)
  mtext(bquote("R"^"2"~"="~.(r_sq)),adj=0.05,line=-4.6)
  mtext(bquote("RMSE = "~.(rmse)),adj=0.05,line=-2)
  mtext(bquote("MAE = "~.(mae)),adj=0.05,line=-3.3)
dev.off() 

# ------------------------------------------------------------------------------
#  Figure S4 PLS Biplots and Coefplots

# load PLS models
pls_models <- list.files("./script/models/", 
                         pattern="pls_", recursive=TRUE,
                         full.names = TRUE)
list_pls_models <- list()
for(pls_model in pls_models){
  temp <- readRDS(pls_model)
  list_pls_models <- append(list_pls_models, list(temp))
}

# plotting
g_pls1_1 <- func_patchwork_plsplot(
  list_pls_models[[1]]$finalModel, 3, 
  list_pls_models[[1]]$finalModel$tuneValue$ncomp, "PLS 1-1")

g_pls1_2 <- func_patchwork_plsplot(
  list_pls_models[[2]]$finalModel, 3, 
  list_pls_models[[2]]$finalModel$tuneValue$ncomp, "PLS 1-2")

g_pls1_3 <- func_patchwork_plsplot(
  list_pls_models[[3]], 3, 
  list_pls_models[[3]]$ncomp, "PLS 1-3", FALSE)

g_pls1_4 <- func_patchwork_plsplot(
  list_pls_models[[4]], 3, 
  list_pls_models[[4]]$ncomp, "PLS 1-4", FALSE)

g_pls2_1 <- func_patchwork_plsplot(
  list_pls_models[[5]]$finalModel, 3, 
  list_pls_models[[5]]$finalModel$tuneValue$ncomp, "PLS 2-1")

g_pls2_2 <- func_patchwork_plsplot(
  list_pls_models[[6]]$finalModel, 3, 
  list_pls_models[[6]]$finalModel$tuneValue$ncomp, "PLS 2-2")

gg_pls <- g_pls1_1 / g_pls1_2 / g_pls1_3 / g_pls1_4 / g_pls2_1/ g_pls2_2   

# export as pdf
ggsave(file="./figures/figure_s4_pls_models.pdf",
       plot=gg_pls,
       device="pdf",
       width=2200, height=4200,
       units="px")

# -----------------------------------------------------------------------------
#  Figure S5; S6 Biplot and standardised coefficients for pls-beta models

plsbeta_models <- list.files("./script/models/", 
                         pattern="plsbeta", recursive=TRUE,
                         full.names = TRUE)
list_plsbeta_models <- list()
for(plsbeta_model in plsbeta_models){
  temp <- readRDS(plsbeta_model)
  list_plsbeta_models <- append(list_plsbeta_models, list(temp))
}

g_plsbeta1_1 <- func_patchwork_plsbetaplot(list_plsbeta_models[[1]], score_scale=2, n_comp=3, "PLS-beta 1-1")
g_plsbeta1_2 <- func_patchwork_plsbetaplot(list_plsbeta_models[[2]], score_scale=2, n_comp=3, "PLS-beta 1-2")
g_plsbeta1_3 <- func_patchwork_plsbetaplot(list_plsbeta_models[[3]], score_scale=2, n_comp=3, "PLS-beta 1-3")
g_plsbeta1_4 <- func_patchwork_plsbetaplot(list_plsbeta_models[[4]], score_scale=2, n_comp=3, "PLS-beta 1-4")
g_plsbeta2_1 <- func_patchwork_plsbetaplot(list_plsbeta_models[[5]], score_scale=7, n_comp=2, "PLS-beta 2-1")
g_plsbeta2_2 <- func_patchwork_plsbetaplot(list_plsbeta_models[[6]], score_scale=7, n_comp=2, "PLS-beta 2-2")

gg_plsbeta_gr <- g_plsbeta1_1 / g_plsbeta1_2 / g_plsbeta1_3 / g_plsbeta1_4 
gg_plsbeta_ng <- g_plsbeta2_1/ g_plsbeta2_2   

gg_plsbeta <- g_plsbeta1_1 / g_plsbeta1_2 / g_plsbeta1_3 / g_plsbeta1_4 /g_plsbeta2_1/ g_plsbeta2_2   


ggsave(file="./figures/figure_s5_plsbeta_models_gr.pdf",
       plot=gg_plsbeta_gr,
       device="pdf",
       width=3200, height=3800,
       units="px")
ggsave(file="./figures/figure_s6_plsbeta_models_ng.pdf",
       plot=gg_plsbeta_ng,
       device="pdf",
       width=2400, height=1900,
       units="px")

# -----------------------------------------------------------------------------
#  Figure S7 RF VarImp plot and Partial dependence plots

# load RandomForest models
rf_models <- list.files("./script/models/", 
                         pattern="rf_", recursive=TRUE,
                         full.names = TRUE)

# variables to process for loop
model_names <- c("RF1-1","RF1-2","RF1-3","RF1-4",
                 "RF2-1","RF2-2")
n_rows <- c(rep(1,5),2)

# drawing graphs for each model
for (i in seq_along(rf_models)){
  temp_rf <- readRDS(rf_models[i])
  save_module_plot(paste0("./figures/rf_plots/", model_names[i], ".png"), 
                   temp_rf, model_names[i], n_rows[i])
}

# stack plots for each model then export
image_files <- list.files("./figures/rf_plots/", pattern=".png", full.names = TRUE)
imgs <- image_read(image_files)
final <- image_append(imgs, stack = TRUE)
image_write(final, "./figures/figure_s7_rf_models.png")
