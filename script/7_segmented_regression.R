# Script for 
#   - Figure 5 time-series plot

library(segmented)
library(tidyverse)
library(jsonlite)
library(ggplot2)
library(ggh4x)
library(patchwork)
library(scales)
library(sf)

gc()
rm(list=ls())
rm(list=ls())

# ------------------------------------------------------------------------------
# load dataset
df_random_sample <- read.csv("./CSV/random_sample_env.csv")

# Load script
source("./script/function/functions_segmented.R")

# Load necessary shapefile
nama_karoo <- st_read("./SHP/NamaKarooDissolved2024.shp") %>%
  st_transform(crs = 4326)
south_africa <- st_read("./SHP/SouthernAfricanCountries.shp") %>%
  filter(NAME=="South Africa")

set.seed(123)  

#-------------------------------------------------------------------------------
# Preprocess df
df_segment <- clean_df_segment(df_random_sample)
# head(df_segment)

# ------------------------------------------------------------------------------
# Run segmented regression
df_segment_nest <- df_segment %>% 
  dplyr::select(vegType,bioregion,random_id,plant_year,bareground) %>% 
  group_by(random_id) %>%
  nest()  # group by pt_id

# Apply 'fit_segmented_model' to rows
# approx 78 seconds for 100 pixels 
# approx 2.85 hours for 15000 pixels

system.time({
  results <- df_segment_nest %>%
    mutate(segmented_model = map(data, fit_segmented_model)) %>%
    unnest(segmented_model)
})

saveRDS(results,"./script/models/segmented_regression_BIC.rds")
# rm(df_segment_nest)

# ------------------------------------------------------------------------------
# Change Point summary
results <- readRDS("./script/models/segmented_regression_BIC.rds")

# N for each vegetation type
df_count <- df_segment %>%
  group_by(vegType) %>%
  distinct(random_id_f) %>% 
  summarise(N = n(), .groups = "drop")

# Extract directions of changes at each of breakpoint
df_change_classification <- classify_change_direction(results)

# summarise the change types by year and vegetation types for plotting
df_change_summary <- summarise_change_pattern(df_change_classification)

# ------------------------------------------------------------------------------
# Figure 5a: Line plot 

# calculate mean per vegtype
df_ts_line_plot <- df_segment %>%
  group_by(bioregion, vegType, plant_year) %>%
  summarise(mean_bareground = mean(bareground, na.rm = TRUE),
            sd_bareground = sd(bareground, na.rm = TRUE),
            .groups = "drop") 

plot_ts_line <- function(dataframe, palette){
  ggplot(dataframe)+
    geom_ribbon(aes(x = plant_year,
                    ymin = mean_bareground-sd_bareground, ymax=mean_bareground+sd_bareground,
                    fill = vegType), alpha=0.25) +
    geom_line(aes(x = plant_year, y = mean_bareground, color = vegType),
              linewidth = 1) +
    scale_y_continuous(labels = function(x) x * 100, limits = c(0,1)) +
    scale_fill_viridis_d(option = palette, begin = 0.3, end = 0.9, drop = TRUE, alpha=0.3) +
    scale_color_viridis_d(option = palette, begin = 0.3, end = 0.9, drop = TRUE) +
    labs(x = "Year", y = "Bare ground (%)") +
    theme_bw() +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5)
    )
  }

ts_line_plot_bl <- df_ts_line_plot %>% 
  filter(bioregion == "Bushmanland") %>% 
  plot_ts_line("turbo")  +
  labs(x = NULL, y = NULL, 
       fill="Bushmanland Bioregion", color="Bushmanland Bioregion", 
       title ="Mean bare ground") 

ts_line_plot_uk <- df_ts_line_plot %>% 
  filter(bioregion == "UpperKaroo") %>% 
  plot_ts_line("plasma")  +
  labs(x = NULL, y = "Bare ground (%)", 
       fill="Upper Karoo Bioregion", color="Upper Karoo Bioregion", 
       title = NULL) 

ts_line_plot_lk <- df_ts_line_plot %>% 
  filter(bioregion == "LowerKaroo") %>% 
  plot_ts_line("mako")+
  labs(x = "Year", y = NULL, 
       fill="Lower Karoo Bioregion", color="Lower Karoo Bioregion", 
       title =NULL) 


ts_line_plot <- ts_line_plot_bl / ts_line_plot_uk / ts_line_plot_lk + plot_layout(ncol=1, guides="collect")

# ------------------------------------------------------------------------------
# Figure 5b: Smoothed slope chart

df_smoothed_slope <- results %>%
    rowwise() %>%
    # open up "fit" column
    mutate(
      slope_df = list({
        fit_vec <- fitted$fit
        year_vec <- fitted$year
        slopes <- c(NA, diff(fit_vec))  # year-to-year differences
        
        tibble(
          random_id = random_id,
          year  = year_vec,
          slope = slopes
        )
      })
    ) %>%
    ungroup() %>%
    dplyr::select(slope_df, vegType, bioregion) %>%
    unnest(slope_df) %>% 
    dplyr::filter(!is.na(slope))
  
# calculate mean per vegtype
df_slope_line_plot <- df_smoothed_slope %>%
  group_by(bioregion, vegType, year) %>%
  summarise(mean_slope = mean(slope, na.rm = TRUE),
            sd_slope = sd(slope, na.rm = TRUE),
            .groups = "drop")   

plot_ts_slope <- function(dataframe, palette){
  ggplot(dataframe)+
    geom_hline(yintercept = 0) +
    geom_ribbon(aes(x = year,
                    ymin = mean_slope-sd_slope, ymax=mean_slope+sd_slope,
                    fill = vegType), alpha=0.25) +
    geom_line(aes(x = year, y = mean_slope, color = vegType),
              linewidth = 1) +
    geom_vline(xintercept = c(1997,2010,2018), linetype = "dotted") +
    scale_y_continuous(labels = function(x) x * 100, limits = c(-0.1,0.1)) +
    scale_fill_viridis_d(option = palette, begin = 0.3, end = 0.9, drop = TRUE, alpha=0.3) +
    scale_color_viridis_d(option = palette, begin = 0.3, end = 0.9, drop = TRUE) +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5)
    )
}

slope_line_plot_bl <- df_slope_line_plot %>%
  filter(bioregion == "Bushmanland") %>% 
  plot_ts_slope("turbo")  +
  labs(x = NULL, y = NULL, 
       fill="Bushmanland Bioregion", color="Bushmanland Bioregion", 
       title ="Slope of change") 

slope_line_plot_uk <- df_slope_line_plot %>% 
  filter(bioregion == "UpperKaroo") %>% 
  plot_ts_slope("plasma")  +
  labs(x = NULL, y = "Annual change in bare ground (%)", 
       fill="Upper Karoo Bioregion", color="Upper Karoo Bioregion", 
       title = NULL) 

slope_line_plot_lk <- df_slope_line_plot %>% 
  filter(bioregion == "LowerKaroo") %>% 
  plot_ts_slope("mako")+
  labs(x = "Year", y = NULL, 
       fill="Lower Karoo Bioregion", color="Lower Karoo Bioregion", 
       title =NULL) 

slope_line_plot <- slope_line_plot_bl + theme(legend.position = "none") +
  slope_line_plot_uk + theme(legend.position = "none") +
  slope_line_plot_lk + theme(legend.position = "none") +
  plot_layout(ncol=1, guides="collect") 

legend <- cowplot::get_legend(slope_line_plot_bl + slope_line_plot_uk + slope_line_plot_lk + plot_layout(ncol=1, guides="collect"))

# ------------------------------------------------------------------------------
# Figure 5c: Break Points chart

df_bioregion_vegtype <- df_ts_line_plot %>% 
  filter(plant_year == 1990) %>% 
  dplyr::select(bioregion, vegType)

df_bioregion_count <- data.frame(
  bioregion = c("Bushmanland", "UpperKaroo", "LowerKaroo"),
  N = c(5000, 5000, 5000)
)

# Format for plotting
df_break_plot <- df_change_summary %>% 
  left_join(df_bioregion_vegtype, by="vegType") %>% 
  filter(change_type =="from_browning_to_greening"|change_type =="from_greening_to_browning") %>% 
  group_by(bioregion,year,change_type) %>% 
  summarise(n = sum(n),
            .groups = "drop") %>% 
  left_join(df_bioregion_count, by="bioregion") %>% 
  mutate(n_ratio = (n/N)*100) %>% 
  mutate(n_ratio=ifelse(change_type =="from_greening_to_browning", -n_ratio, n_ratio)) 

breakpoint_plot <- df_break_plot %>% 
  ggplot(aes(x = year, y = n_ratio)) +
  geom_col(data = filter(df_break_plot, change_type =="from_greening_to_browning"),fill="forestgreen") +
  geom_col(data = filter(df_break_plot, change_type =="from_browning_to_greening"),fill="saddlebrown", position = position_stack(reverse = TRUE)) +
  scale_y_continuous(labels = function(x) abs(x), limits=c(-30,30))+
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = c(1997,2010,2018), linetype = "dotted") +
  annotate("text", label = "Greening Peak" , x = 1998, y = -28, size = 3, colour = "grey20")+
  annotate("text", label = "Browning Peak", x = 1998, y = 28, size = 3, colour = "grey20")+
  theme_bw()+
  labs(x = "Year", y = "Breakpoint occurrence (%)",
       title = "Breakpoints counts")+
  facet_wrap(~bioregion, ncol=1) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

gg_timeseries <- (ts_line_plot) | (slope_line_plot) | (breakpoint_plot)|legend

ggsave("./figures/figure5_time-series_analysis.pdf",
       plot = gg_timeseries,
       width = 3600, height = 1800, units="px")
