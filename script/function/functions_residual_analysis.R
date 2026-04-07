#-------------------------------------------------------------------------------
# Functions  
#-------------------------------------------------------------------------------

# Preparing validation data - same as validation.R
func_pred_valid <- function(model, dataframe){
  output <- dataframe %>% 
    dplyr::select(BareGround, Biome, BioRegion, Year, Month) %>% 
    rownames_to_column(var="ID") %>%
    mutate(ID=as.character(ID))
  
    if(class(model)=="train"){                         # caret
      temp_df <- dataframe %>% 
        dplyr::select(all_of(c("BareGround",temp_model$finalModel$xNames))) %>% 
        na.omit()
      temp_pred <- predict(model, temp_df)
    } else if (class(model)=="plsRbetamodel") {        # plsR-beta
      temp_df <- dataframe %>% 
        dplyr::select(all_of(c("BareGround",colnames(model$dataX)))) %>% 
        na.omit()
      temp_pred <- predict_plsRbeta(model, temp_df)*100
    } else {
      temp_df <- dataframe %>% 
        dplyr::select(all_of(
          colnames(model$model)
        )) %>% 
        na.omit()
      temp_pred <- predict(model, temp_df)[,,2]  # mvr or pls
    }
    temp_colnames <- c(colnames(output), "predicted")
    temp_pred_df <- data.frame(temp = temp_pred) %>% 
      mutate(ID=row.names(temp_df))
    output <- output %>%
      left_join(temp_pred_df, by="ID")
    colnames(output) <- temp_colnames

  output <- output %>%
    dplyr::select(-ID)
   return(output)
}


# ------------------------------------------------------------------------------
# Function to plot standardized residuals vs predicted values
# Inputs:
#   - df: A data frame containing at least the following columns:
#         predicted, std_residuals, model, and one of the color_column
#   - color_column: A string specifying the column to color points by (e.g., "Month", "Year")
#   - color_option: A string specifying the viridis color palette option (e.g., "D", "H")
#   - color_label: A string specifying the label of the color legend
# Output:
#   - A ggplot object
# ------------------------------------------------------------------------------
plot_residuals_by_covariate <- function(df, color_column, color_option = "D", color_label = NULL, color_lim) {
  # color_label <- color_label %||% color_column  # fallback if color_label is NULL
  
  ggplot(df, aes(x = predicted, y = residuals,
                 color = .data[[color_column]],
                 shape = data_type)) +  
    geom_point(size = 3, alpha = 0.25) +
    scale_color_viridis_c(option = color_option, name = color_label,
                          limits = color_lim, oob = scales::squish) +
    scale_shape_manual(
      values = c("photo_train" = 16, "photo_valid" = 17),  
      name = "Dataset",
      labels=c("Training", "Validation")
    ) +
    geom_hline(yintercept = 0, color = "black", show.legend = FALSE) +
    ylim(c(-60, 60)) +
    xlim(c(0, 100)) +
    labs(x = "Predicted bare ground (%)", y = "Differences between \n observed and predicted values (%)", color = color_label) +
    # facet_wrap(~ model, nrow = 3) +
    theme_bw()
}


plot_residuals_by_biome <- function(df) {
  df %>% mutate(Biome_BioRegion = case_when(
    Biome %in% c("Albany Thicket", "Azonal Vegetation", "Desert", "Grassland") ~ Biome,
    TRUE ~ BioRegion
  )) %>%
    mutate(colour = color_mapping[Biome_BioRegion]) %>%
    # plotting
    ggplot(aes(x = predicted, y = residuals, color = colour, shape = data_type)) +
    geom_point(size = 3, alpha = 0.25) +
    scale_shape_manual(
      values = c("photo_train" = 16, "photo_valid" = 17),  
      name = "Dataset",
      labels = c("Training", "Validation")
    ) +
    scale_color_identity(
      guide = "legend",
      name = "Biome/Bioregion",
      breaks = names(label_mapping),
      labels = label_mapping
    ) +
    geom_hline(yintercept = 0, color = "black", show.legend = FALSE) +
    ylim(-60, 60) +
    xlim(0, 100) +
    labs(x = "Predicted bare ground (%)", y = "Differences between \n observed and predicted values (%)", color = "Year") +
    theme_bw()
}

