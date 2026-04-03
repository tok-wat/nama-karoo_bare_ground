# Script for 
#   - Figure 6  map and ridge plot
#   - Figure S4 Composite annual map

library(sf)
library(tidyverse)
# library(raster)
library(terra)
library(tidyterra)
library(ggplot2)
library(ggspatial)
library(ggridges)
library(RColorBrewer)
library(patchwork)


##################################
rm(list=ls())

gc(reset = TRUE)
gc(reset = TRUE)
##################################
# Load a tiff
pred_annual_200 <- rast("./TIFF/pred_annual_200m.tif")
pred_annual_200[pred_annual_200 == 0] <- NA

# Load necessary shapefile
nama_karoo <- st_read("./SHP/NamaKarooDissolved2024.shp") %>%
  st_transform(crs = 4326)
south_africa <- st_read("./SHP/SouthernAfricanCountries.shp") %>%
  filter(NAME=="South Africa")

#-------------------------------------------------------------------------------
# Figure 6A: Maps

bare_1990 <- pred_annual_200[[1]]
bare_1997 <- pred_annual_200[[8]]  # 1997 - 1990 + 1 = 8
bare_2010 <- pred_annual_200[[21]]
bare_2018 <- pred_annual_200[[29]]
bare_2024 <- pred_annual_200[[35]]

plot_map <- function(spat_rast, year){
  ggplot() +
  # geom_sf(data = nama_karoo, fill = "grey40", color = NA, alpha = 0.3) +
  geom_spatraster(data = spat_rast) +
  geom_sf(data = south_africa, fill = NA) +
  coord_sf(xlim = c(18, 28), ylim = c(-34, -26)) +
  annotate("text", x = 26.5, y = -26.3, label = year, size = 4) +
  scale_fill_gradientn(
      labels = function(x) x / 100,
      colors = rev(brewer.pal(11, "BrBG")),
      limits = c(0, 10000),
      oob = scales::squish,
      na.value = NA,
      name = "Bare Ground (%)"
    )+
  theme_bw()+
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
}

map_1990 <- plot_map(bare_1990, 1990)
map_1997 <- plot_map(bare_1997, 1997)
map_2010 <- plot_map(bare_2010, 2010)
map_2018 <- plot_map(bare_2018, 2018)
map_2024 <- plot_map(bare_2024, 2024)

maps <- map_1990 + map_1997 + map_2010 + map_2018 + map_2024 +
  plot_layout(nrow = 2, guides="collect")

# ------------------------------------------------------------------------------
# Figure 6B: Ridge Plot

df_segment <- readRDS("./script/models/segmented_regression_BIC.rds")

df_ridge <- df_segment %>% 
  dplyr::select(random_id, data) %>% 
  unnest(data) %>% 
  filter(plant_year %in% c(1990, 1997, 2010, 2018, 2024)) %>% 
  group_by(random_id) %>%
  mutate(
    bareground = bareground*100,
    bareground_diff = bareground - bareground[plant_year == 1990]) %>% 
  ungroup() %>% 
  filter(plant_year != 1990)

# ridge plot for cover data
plot_ridge <-function(dataframe, palette){
  ggplot(dataframe, 
         aes(x = bareground_diff, 
             y = factor(plant_year, levels = c("2024", "2018","2010","1997")), 
             fill = vegType,
             color = vegType)) +
    geom_density_ridges(
      scale = 0.8,
      rel_min_height = 0.02,
      bandwidth = bw_value,
      position="identity",
      linewidth = 1
    ) +
    scale_fill_viridis_d(option = palette, begin = 0.3, end = 0.9, drop = TRUE, alpha=0.3) +
    scale_color_viridis_d(option = palette, begin = 0.3, end = 0.9, drop = TRUE) +
    coord_cartesian(xlim = c(-40, 40)) +
    scale_y_discrete(expand = expansion(mult = c(0, 0.3))) + # space between plots and label
    geom_vline(xintercept = 0, linetype = "dotted") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
}
bw_value <- 1.7

ridge_bushmanland <- df_ridge %>% 
  filter(bioregion == "Bushmanland") %>% 
  plot_ridge("turbo") +
  labs(x = NULL, y = "Density", 
       fill="Bushmanland Bioregion", color="Bushmanland Bioregion", 
       title ="Bushmanland") 

ridge_upperkaroo <- df_ridge %>% 
  filter(bioregion == "UpperKaroo") %>% 
  plot_ridge("plasma") +
  labs(x = "Change in bare ground from 1990 (%)", y = NULL, 
       fill="Upper Karoo Bioregion", color = "Upper Karoo Bioregion",
       title = "Upper Karoo")

ridge_lowerkaroo <- df_ridge %>% 
  filter(bioregion == "LowerKaroo") %>% 
  plot_ridge("mako") +
  labs(x=NULL, y=NULL, 
       fill="Lower Karoo Bioregion", color="Lower Karoo Bioregion",
       title = "Lower Karoo")

ridges <- ridge_bushmanland + ridge_upperkaroo + ridge_lowerkaroo +
   plot_layout(guides="collect") 


# Export
ggsave(file="./figures/figure6a_maps.pdf",
       plot=maps,
       device="pdf",
       width=2800, height=1600,
       units="px")

ggsave(file="./figures/figure6b_ridge_plot.pdf",
       plot=ridges,
       device="pdf",
       width=2800, height=1000,
       units="px")

# ------------------------------------------------------------------------------
# Figure S9: Composite annual maps

# Create labels for each facet (plant year) to display as text
label_df <- data.frame(
  lyr = names(pred_annual_200),
  label = 1990:2024,
  x = 25.3,   # x-coordinate for the text
  y = -27.2    # y-coordinate for the text
)

# Create the plot
p <- ggplot() +
  #  geom_sf(data = nama_karoo, fill = "grey40", color = NA, alpha = 0.6) +
  geom_spatraster(data = pred_annual_200) +
  geom_sf(data = south_africa, fill = NA) +
  coord_sf(xlim = c(18, 27), ylim = c(-34, -26.5)) +
  scale_fill_gradientn(
    labels = function(x) x / 100,
    colors = rev(brewer.pal(11, "BrBG")),
    limits = c(0, 10000),
    oob = scales::squish,
    na.value = NA,
    name = "Bare Ground (%)"
  )+
  facet_wrap(~ lyr, ncol = 5) +
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 5,            # Text size (adjust as needed)
    color = "black"      # Text color (adjust as needed)
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.text = element_blank(),         # Hide facet strip text (labels)
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(fill = "Fraction of bareground (%)")

# Export
ggsave("./figures/figure_s9_annual_bareground_composite_map.png", p,
       width=10, height=16)