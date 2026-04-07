# ------------------------------------------------------------------------------
# Setting for plotting

# labels and colour
vi_labels <- c(
  exg = "ExG",
  swirRatio = "SWIR ratio",
  tc_bright = "TC Brightness",
  tc_green = "TC Greeness",
  tc_wet = "TC Wetness",
  elevation = "Elevation",
  potassium = "Potassium",
  soil_depth = "Soil depth",
  stone = "Stone"
)

vi_colours <- c(
  ndmi = "blue3",
  tc_wet = "blue3",
  bsi = "red3",
  tc_bright = "red3",
  swirRatio = "orange",
  satvi = '#808000',
  elevation = "grey",
  potassium = "grey",
  soil_depth = "grey",
  stone = "grey"
)

# ------------------------------------------------------------------------------
# Function to compute the proportion of X-variance explained by each component
# Input
# Output

func_exp_var <- function(plsRbeta_object) {
  X <- plsRbeta_object$dataX
  T_scores <- plsRbeta_object$tt
  P_loadings <- plsRbeta_object$pp
  
  # Total sum of squares of X
  SSX_total <- sum(scale(X)^2)
  
  # Sum of squares explained by each component
  SSXcomp <- sapply(1:ncol(T_scores), function(k) {
    t_k <- T_scores[, k]
    p_k <- P_loadings[, k]
    X_k_hat <- t_k %*% t(p_k)  # reconstructed X (n × p)
    sum(X_k_hat^2)            # variance explained by component k
  })
  
  # Proportion of explained variance
  explained_X <- SSXcomp / SSX_total
  
  # Return as a data frame
  df_x_var <- data.frame(
    Comp = paste0("Comp_", 1:length(explained_X)),
    Expl_Var_X = explained_X,
    Cumulative = cumsum(explained_X)
  )
  return(df_x_var)
}

# ------------------------------------------------------------------------------
# Function to plot PLS-beta regression biplot with score points and variable arrows
plot_biplot <- function(plsRbeta_object, 
                        label_scale = 1.2, score_scale, components = 1:2,
                        show_legend=TRUE) {
  # Extract and format scores
  S <- plsRbeta_object$tt %>%
    as.data.frame() %>%
    rename_with(~ str_replace_all(., " ", "")) %>%
    mutate(.outcome = plsRbeta_object$dataY * 100)
  
  # Extract pseudo R2
  pseudoR2 <- plsRbeta_object$pseudo.R2[1,1]
  pseudoR2 <- c(pseudoR2, plsRbeta_object$pseudo.R2%>% as.vector() %>% diff())
  
  # Extract and format loadings, compute label coordinates
  P_loadings <- plsRbeta_object$pp %>%
    as.data.frame() %>%
    rename_with(~ str_replace_all(., " ", "")) %>%
    mutate(across(starts_with("Comp_"), ~ .x * score_scale, .names = "{.col}")) %>%
    mutate(across(starts_with("Comp_"), ~ .x * label_scale, .names = "label_{.col}")) %>%
    rownames_to_column(var = "variable") %>% 
    mutate(variable = str_remove(variable, "_boxcox"))  
  
  # Get explained variance for axis labels
  # df_x_var <- func_exp_var(plsRbeta_object)
  # xlab_str <- bquote(
  #   "Comp " * .(components[1]) *
  #     " (R"[x]^2 * " = " * .(sprintf("%.2f", round(df_x_var[components[1], 2], 3))) * ", " *
  #     "pseudo R"[y]^2 * " = " * .(sprintf("%.2f", round(pseudoR2[components[1]] , 3))) *")"
  # )
  # 
  # ylab_str <- bquote(
  #   "Comp " * .(components[2]) *
  #     " (R"[x]^2 * " = " * .(sprintf("%.2f", round(df_x_var[components[2], 2] , 3))) * ", " *
  #     "pseudo R"[y]^2 * " = " * .(sprintf("%.2f", round(pseudoR2[components[2]], 3))) * ")"
  # )
  xlab_str <- paste("Comp", components[1])
  ylab_str <- paste("Comp", components[2])
  
  # Specify column names used for plotting
  xcomp <- colnames(P_loadings)[components[1] + 1]
  ycomp <- colnames(P_loadings)[components[2] + 1]
  xcomp_vi_label <- paste0("label_", xcomp)
  ycomp_vi_label <- paste0("label_", ycomp)
  
  # Determine axis limits based on scores, arrows, and labels
  axis_lim <- c(
    x = max(abs(c(S[[xcomp]], P_loadings[[xcomp]], P_loadings[[xcomp_vi_label]]))),
    y = max(abs(c(S[[ycomp]], P_loadings[[ycomp]], P_loadings[[ycomp_vi_label]])))
  ) * 1.05
  
  # Construct the biplot
  p <- ggplot() +
    # Draw axes lines
    geom_hline(yintercept = 0, color = "gray80") +
    geom_vline(xintercept = 0, color = "gray80") +
    
    # Plot score points
    geom_point(data = S, aes(x = !!sym(xcomp), y = !!sym(ycomp), color = .outcome),
               size = 2, alpha = 0.6,
               show.legend=show_legend) +
    
    # Add variable arrows
    geom_segment(data = P_loadings,
                 aes(x = 0, y = 0, xend = !!sym(xcomp), yend = !!sym(ycomp)),
                 arrow = arrow(length = unit(0.2, "cm")),
                 linewidth = 1.5, alpha = 0.5) +
    
    # Add labels for variables
    geom_text(data = P_loadings,
              aes(x = !!sym(xcomp_vi_label), y = !!sym(ycomp_vi_label), label = variable),
              size = 4) +
    
    # Set symmetric axis limits and secondary axes for variable loadings
    scale_x_continuous(
      limits = c(-axis_lim["x"], axis_lim["x"]),
      sec.axis = sec_axis(~ . / score_scale,
                          name = paste0("Variable loading (Comp ", components[1], ")"))
    ) +
    scale_y_continuous(
      limits = c(-axis_lim["y"], axis_lim["y"]),
      sec.axis = sec_axis(~ . / score_scale,
                          name = paste0("Variable loading (Comp ", components[2], ")"))
    ) +
    
    # Theme and labels
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)) +
    labs(x = xlab_str, y = ylab_str) +
    
    # Gradient color scale for outcome variable
    scale_color_gradientn(
      colors = rev(brewer.pal(11, "BrBG")),
      limits = c(0, 100),
      name = "Bare Ground (%)"
    )
  
  return(p)
}

# combine biplot and standardised coefficients
func_patchwork_plsbetaplot <- function(plsbeta_obj, score_scale=2, n_comp=2, title){

  temp_coef_plot <- process_std_coef(plsbeta_obj$Std.Coeffs, vi_labels, vi_colours) %>% 
    plot_std_coef()  
  
  if(n_comp == 2){
    temp_biplot <- plot_biplot(plsbeta_obj, label_scale=1.2, 
                             score_scale, components = 1:2) +
    labs(title=title)
    g <- temp_biplot +  temp_coef_plot + plot_layout(ncol = 2, widths = c(1, 1))
  }else if(n_comp == 3){
    temp_biplot_1_2 <- plot_biplot(plsbeta_obj, label_scale=1.2, 
                               score_scale, components = 1:2, show_legend = FALSE) +
      labs(title=title)  
    temp_biplot_1_3 <- plot_biplot(plsbeta_obj, label_scale=1.2, 
                               score_scale, components = c(1,3)) 
    g <- temp_biplot_1_2 +  temp_biplot_1_3 + temp_coef_plot + 
      plot_layout(ncol = 3, widths = c(1, 1, 1))
  }
  return(g)
}
# ------------------------------------------------------------------------------
# Function to draw biplot from a PLS-beta regression model

plot_pls_biplot <- function(pls_obj, arrows_scale=3, title="PLS Model"){
  # --- extract X scores ---
  x_scores <- pls_obj$scores
  scores_plot <- as.data.frame(x_scores[, 1:2]) %>%
    rename(Comp1 = 1, Comp2 = 2)
  
  # --- calculate correlation with X ---
  X_mat <- as.matrix(pls_obj$model[, -1])
  cor_X_df <- cor(X_mat, x_scores) %>%
    as.data.frame() %>%
    rownames_to_column(var = "Variable") %>%
    mutate(Type = "X",
           Variable = str_remove(Variable, "_boxcox"))  
  
  
  # --- calculate correlation with Y ---
  Y_vec <- pls_obj$model[, 1]
  cor_Y_df <- cor(Y_vec, x_scores) %>%
    as.data.frame() %>%
    mutate(Variable ="Y",
           Type = "Y")
  
  # --- combine X and Y ---
  cor_all_df <- bind_rows(cor_X_df, cor_Y_df)
  
  # --- arrows ---
  cor_plot_df <- cor_all_df %>%
    rename(Comp1 = 2, Comp2 = 3) %>%  
    mutate(
      xend = Comp1 * arrows_scale,
      yend = Comp2 * arrows_scale
    )
  
  # --- setting drawing ranges ---
  xmax <- ceiling(max(scores_plot$Comp1))
  xmin <- floor(min(scores_plot$Comp1))

  # --- ggplot ---
  p <- ggplot() +
    geom_point(data = scores_plot, aes(x = Comp1, y = Comp2), color = "grey20", alpha = 0.5) +
    geom_segment(data = cor_plot_df, aes(x = 0, y = 0, xend = xend, yend = yend, color = Type),
                 arrow = arrow(length = unit(0.2, "cm")), size = 1) +
    geom_text(data = cor_plot_df, aes(x = xend * 1.2, y = yend * 1.2, label = Variable, color = Type),
              size = 4) +
    scale_color_manual(values = c("X" = "red", "Y" = "blue")) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "Comp 1", y = "Comp 2", title=title) +
    xlim(c(xmin, xmax)) 
  
  return(p)
}
plot_pls_vip <- function(pls_obj, n_comp){
  g <- VIP(pls_obj, n_comp) %>% 
    as.data.frame() %>%
    rownames_to_column("vi") %>%
    rename(VIP = ".") %>%
    mutate(
      vi = str_remove(vi, "_boxcox"),  
      colour = coalesce(vi_colours[vi], "green4"),
      vi = coalesce(vi_labels[vi],  toupper(vi)),
      vi = fct_reorder(vi, VIP)) %>% 
    ggplot()+
    geom_col(aes(x=vi, y=VIP, fill=colour))+
    scale_fill_identity() +
    coord_flip()+
    ylim(c(0,1.5))+ 
    xlab("Vegetation Indices") +
    theme_bw() +
    theme(axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8))
  return(g)
}
plot_pls_std_coef <- function(pls_obj, scale_X = TRUE){
  reg_beta <- coef(pls_obj) %>% 
    as.vector()
  
  X <- pls_obj$model[,-1] # predictors
  Y <- pls_obj$model[,1]  # respons variable
  
  # get Std dev then standardise
  sd_Y <- sd(Y)
  if(isTRUE(scale_X)){
    sd_X <- apply(X, 2, sd)
    std_beta <- reg_beta * sd_X / sd_Y
  }else{
    std_beta <- reg_beta / sd_Y
    names(std_beta) <- colnames(X)
  }
  
  g <- std_beta %>% 
    as.data.frame() %>%
    rownames_to_column("vi") %>%
    rename(std_coef = ".") %>%
    mutate(
      vi = str_remove(vi, "_boxcox"),
      colour = coalesce(vi_colours[vi], "green4"),
      vi = coalesce(vi_labels[vi],  toupper(vi)),
      vi = fct_reorder(vi, std_coef)
    ) %>% 
    plot_std_coef() +
    theme(axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8))
  return(g)
}

func_patchwork_plsplot <- function(pls_obj, arrows_scale=3, n_comp, title, scale_X=TRUE){
  temp_biplot    <- plot_pls_biplot(pls_obj, arrows_scale, title)
  temp_vip_plot  <- plot_pls_vip(pls_obj, n_comp)
  temp_coef_plot <- plot_pls_std_coef(pls_obj, scale_X)
  
  g <- temp_biplot + temp_vip_plot + temp_coef_plot + plot_layout(ncol=3, width=c(1.2,1,1))
  return(g)
}

# ------------------------------------------------------------------------------
# Function to plot standardised coefficients from a PLS-beta regression model


# processing function
process_std_coef <- function(std_coef_mat, vi_labels, vi_colours) {
  std_coef_mat[-1, 1] %>%
    as.data.frame() %>%
    rownames_to_column("vi") %>%
    rename(std_coef = ".") %>%
    mutate(
      vi = str_remove(vi, "_boxcox"),
      colour = coalesce(vi_colours[vi], "green4"),
      vi = coalesce(vi_labels[vi],  toupper(vi)),
      vi = fct_reorder(vi, std_coef)
    )
}

# ggplot function
plot_std_coef <- function(df) {
  ggplot(data = df) +
    geom_hline(yintercept = seq_len(nrow(df)), linetype = "dotted", alpha = 0.5) +
    geom_vline(xintercept = 0, alpha = 0.5) +
    geom_segment(aes(x = std_coef, xend = 0, y = vi, yend = vi, color = colour),
                 linewidth = 1.5, alpha = 0.7) +
    geom_point(aes(x = std_coef, y = vi, color = colour),
               shape = 16, size = 2.5, stroke = 2) +
    scale_color_identity() +
    xlab("Standardised coefficients") +
    ylab("Vegetation Indices") +
    xlim(-1.0, 1.0) +
    theme_bw()
}

# ------------------------------------------------------------------------------
# Function to plot variable importance and partial dependence from a RandomForest model

# plotting variable importance
plot_var_imp <- function(rf_obj, title){

  # -- extract %IncMSE
  df_imp_plot <- importance(rf_obj, type=1) %>% 
    as.data.frame() %>% 
    rownames_to_column(var="vi") %>% 
    arrange(`%IncMSE`) %>% 
    mutate(
      colour = coalesce(vi_colours[vi], "green4"),
      vi = coalesce(vi_labels[vi],  toupper(vi)),
      vi = fct_reorder(vi, `%IncMSE`)
    )
  
  par(mar = c(4, 8,  2, 2),
      mgp = c(6, 1, 0))
  
  bp <- barplot(`%IncMSE` ~ vi, data=df_imp_plot,   # extract bars positions
                plot = FALSE)
  # -- empty plot --
  plot(NA, xlim = c(0, 70),
       ylim = range(bp)+c(-0.5, 0.5),
       xlab = "", ylab = "Predictors",
       yaxt = "n", bty = "n")
  mtext(side=3, text = title, line=0, adj=0)
  mtext(side=1, text = "%IncMSE", line=2.5, adj=0.5, cex=0.8)
  
  # -- draw plot -- 
  abline(v = seq(10, 60, 10), col = "grey80", lty = 1)
  barplot(`%IncMSE` ~ vi, data=df_imp_plot,
          col = df_imp_plot$colour,
          border = NA,
          las = 1, horiz = TRUE,
          xlim = c(0, 50),
          add = TRUE)
  box(lty = 1)
}

# partial dependence plot
plot_pdp <- function(caret_rf_obj){
  # -- margins
  par(mar = c(4, 3.5, 2, 1),
      mgp = c(2, 0.8, 0))
  
  # -- extract variable importance to order plots 
  temp_df <- importance(caret_rf_obj$finalModel, type=1) %>% 
    as.data.frame() %>% 
    rownames_to_column(var="vi") %>% 
    arrange(desc(`%IncMSE`))
  
  # -- setting to run for loop
  temp_pred_vars <- temp_df$vi
  labels <- data.frame(vi=temp_pred_vars) %>%
    mutate(vi = coalesce(vi_labels[vi],  toupper(vi))) %>% 
    unlist()
  
  # -- use do.call for plotting for each predictor (normal for loop didn't work)
  for(i in seq_along(temp_pred_vars)){
    var <- temp_pred_vars[i]
    lbl <- labels[i]
    do.call(partialPlot, list(
      x = caret_rf_obj$finalModel,
      pred.data = caret_rf_obj$trainingData,
      x.var = as.name(var),  
      rug = TRUE,
      col = "blue",
      xlab = lbl,
      ylab = "Partial dependence",
      main = ""
    ))
  }
}

# For each model, draw var imp plot and partial dependence plot, combine them,
# then export as a png file.
save_module_plot <- function(file_name, caret_rf_obj, title, n_row=1){
  # layout
  if(n_row == 2){
    module <- matrix(c(
      1, 2, 3, 4, 5,
      1, 6, 7, 8, 9 
    ), nrow = 2, byrow = TRUE)
    height <- 600
  }else if(n_row == 1){
    module <- matrix(c(
      1, 2, 3, 4, 5
    ), nrow = 1, byrow = TRUE)
    height <- 300
  }
  
  png(file_name, width=2000, height=height, res=200)  
  
  # left : right ：
  layout(module, widths = c(3, 1, 1, 1, 1)*1, heights = c(1,1))
  
  # var imp plot 
  par(mar=c(5,4,2,1))  
  plot_var_imp(caret_rf_obj$finalModel, title)
  
  # partial dependence plots
  plot_pdp(caret_rf_obj)
  
  dev.off()
}


# ------------------------------------------------------------------------------
# Helper functions to prepare maps for Figure 4

# Function to rescaling the coordinates 
rescale_coords <- function(df) {
  df %>%
    mutate(
      x = 100 * (x - min(x)) / (max(x) - min(x)),
      y = 100 * (y - min(y)) / (max(y) - min(y))
    )
}

# Prepare geo-tiff for visualisation
read_and_crop_raster <- function(raster_path, shapefile_sf) {
  r <- raster(raster_path)
  r_proj <- projectRaster(r, crs = CRS("+init=EPSG:32734"), method = "bilinear")
  # cropping raster with shapefile
  cropped <- mask(crop(r_proj, extent(shapefile_sf)), as(shapefile_sf, "Spatial"))
  df <- as.data.frame(cropped, xy = TRUE)
  colnames(df) <- c("x", "y", "value")
  # rescale
  rescale_coords(df)
}

# Prepare single band RGB geo-tiff for visualisation
prepare_rgb_raster <- function(rgb_path, shapefile_sf) {
  rgb_stack <- stack(rgb_path)
  cropped <- mask(crop(rgb_stack, extent(shapefile_sf)), as(shapefile_sf, "Spatial"))
  df <- as.data.frame(cropped, xy = TRUE)
  colnames(df)[3:5] <- c("R", "G", "B")
  df <- df %>%
    filter(!is.na(R) & !is.na(G) & !is.na(B)) %>%
    mutate(hex = rgb(R, G, B, maxColorValue = 255)) %>%
    dplyr::select(x, y, hex)
  rescale_coords(df)
}

# ------------------------------------------------------------------------------
# ggplot functions for draw sub-plot of Figure 4

# plotting Google Earth Image
plot_ge <- function(df1, annotate){
  p <- ggplot(df1) +
    geom_raster(aes(x = x, y = y, fill = hex)) +
    scale_fill_identity() +
    coord_fixed(xlim = c(0, 100), ylim = c(0, 100), expand = FALSE) +
    theme_void() +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
  annotate("text", x = 90, y = 5, label = annotate, size = 6)
  return(p)
}
# plotting the bare ground fraction rasters
plot_rast <- function(df2, annotate){
  p <- ggplot(df2) +
    geom_raster(aes(x = x, y = y, fill = value)) +
    #scale_fill_viridis_c(option = "magma", limits = c(0, 100), oob = scales::squish,, direction = -1) +
    scale_fill_gradientn(
      colors = rev(brewer.pal(11, "BrBG")),
      limits = c(0, 100),
      oob = scales::squish,
      name = "Bare Ground (%)"
    )+
    coord_fixed(xlim = c(0, 100), ylim = c(0, 100), expand = FALSE) +
    theme_void() +
    theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
    annotate("text", x = 90, y = 5, label = annotate, size = 6)
  return(p)
}
