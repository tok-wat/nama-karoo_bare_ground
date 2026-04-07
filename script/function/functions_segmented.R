# preprocessing df_random_sample for this part of the analysis

clean_df_segment <- function(df){
  cleaned <- df %>%
    dplyr::select(starts_with("bareground"),bioregion, vegType, .geo) %>% 
    mutate(
      # extract coordinates
      parsed = map(.geo, ~ fromJSON(.)$coordinates),
      lon    = map_dbl(parsed, 1),
      lat    = map_dbl(parsed, 2)
    ) %>%
    dplyr::select(-c(.geo, parsed)) %>%
    rowid_to_column("random_id") %>%
    na.omit() %>%
    pivot_longer(cols=starts_with("bareground"),
                 names_to = "plant_year",
                 values_to = "bareground") %>% 
    separate(plant_year, into = c("prefix", "plant_year", "suffix"), sep = "_") %>% 
    mutate(random_py_id=paste0(random_id,"_",plant_year)) %>%
    # convert to 0 to 1 for beta regression
    mutate(bareground=scales::rescale(bareground, to=c(0.0001, 0.9999), from=range(0,10000))) %>% 
    dplyr::select(-c(prefix, suffix)) %>% 
    mutate(
      # labeling bioregion 
      bioregion = recode_factor(bioregion,
                                `1` = "Bushmanland",
                                `2` = "UpperKaroo",
                                `3` = "LowerKaroo"
      ),
      # labeling vegetation type 
      vegType = recode_factor(vegType,
                              `1`  = "Western Upper Karoo",
                              `2`  = "Upper Karoo Hardeveld",
                              `3`  = "Northern Upper Karoo",
                              `4`  = "Eastern Upper Karoo",
                              `5`  = "Gamka Karoo",
                              `6`  = "Eastern Lower Karoo",
                              `7`  = "Albany Broken Veld",
                              `10` = "Lower Gariep Broken Veld",
                              `11` = "Blouputs Karroid Thornveld",
                              `12` = "Bushmanland Arid Grassland",
                              `13` = "Bushmanland Sandy Grassland",
                              `14` = "Kalahari Karroid Shrubland",
                              `15` = "Bushmanland Basin Shrubland"
      ),
      # ordering factors
      vegType   = factor(vegType, levels = c(
        "Lower Gariep Broken Veld","Blouputs Karroid Thornveld","Bushmanland Arid Grassland",
        "Bushmanland Sandy Grassland","Kalahari Karroid Shrubland","Bushmanland Basin Shrubland",
        "Western Upper Karoo","Upper Karoo Hardeveld","Northern Upper Karoo", "Eastern Upper Karoo",
        "Gamka Karoo", "Eastern Lower Karoo", "Albany Broken Veld"
      )),
      bioregion = factor(bioregion, levels = c("Bushmanland", "UpperKaroo", "LowerKaroo")),
      # factorise plant_year and id
      plant_year = as.numeric(plant_year),
      plant_year_f = factor(plant_year),
      random_id_f  = factor(random_id)
    ) %>%
    group_by(bioregion) %>%
    mutate(random_id = as.integer(random_id)) %>%
    filter(random_id %in% sample(unique(random_id), 5000)) %>%  # select 5000 pts per bioregion
    ungroup()
  
  return(cleaned)
}

# function to fit segmented regression with k =1 to 10, and choose the model 
# with the lowest AIC without error in the process
fit_segmented_model <- function(data) {
  
  lm_fit <- lm(bareground ~ plant_year, data = data)
  
  best_model <- lm_fit
  breakpoints <- list(NA)
  # best_aic <- AIC(lm_fit)
  best_bic <- BIC(lm_fit)
  best_k <- 0
  fitted <- list(fit=fitted(lm_fit), year=data$plant_year)
  
  for (k in 1:10) {
    warning_flag <- FALSE
    seg_fit <- NULL
    
    #  catch errors and warning
    tryCatch({
      # flag warning 
      withCallingHandlers({
        seg_fit <- segmented(lm_fit, seg.Z = ~plant_year, npsi = k, control = seg.control(display = FALSE))
      }, warning = function(w) {
        warning_flag <<- TRUE
        invokeRestart("muffleWarning")  # do not display warnings
      })
    }, error = function(e) {
      # ignore errors
      warning_flag <<- TRUE
    })
    
  #   if (!warning_flag && !is.null(seg_fit)) {
  #     model_aic <- AIC(seg_fit)
  #     if (model_aic < best_aic) {
  #       best_model <- seg_fit
  #       breakpoints <- seg_fit$psi[,"Est."]
  #       best_aic <- model_aic
  #       best_k <- k
  #       fitted <- list(fit=as.vector(best_model$fitted.values),
  #                      year=data$plant_year)
  #     }
  #   }
  # }
  # best_fit <- tibble(vegType=data$vegType[1],
  #                    bioregion=data$bioregion[1],
  #                    breakpoints = list(breakpoints), 
  #                    AIC = best_aic, 
  #                    k = best_k,
  #                    fitted = list(fitted))
      if (!warning_flag && !is.null(seg_fit)) {
        model_bic <- BIC(seg_fit)
        if (model_bic < best_bic) {
          best_model <- seg_fit
          breakpoints <- seg_fit$psi[,"Est."]
          best_bic <- model_bic
          best_k <- k
          fitted <- list(fit=as.vector(best_model$fitted.values),
                         year=data$plant_year)
        }
      }
    }
    best_fit <- tibble(vegType=data$vegType[1],
                       bioregion=data$bioregion[1],
                       breakpoints = list(breakpoints),
                       BIC = best_bic,
                       k = best_k,
                       fitted = list(fitted))
  return(best_fit)
}

# # process the output of segmented regressions and summarise it to a dataframe
# summarise_slope_by_year <- function(tibble_results){
#   df_processed <- tibble_results %>%
#     rowwise() %>%
#     # open up "fit" column
#     mutate(
#       slope_df = list({
#         fit_vec <- fitted$fit
#         year_vec <- fitted$year
#         slopes <- c(NA, diff(fit_vec))  # year-to-year differences
#         
#         tibble(
#           random_id = random_id,
#           year  = year_vec,
#           slope = slopes
#         )
#       })
#     ) %>%
#     ungroup() %>%
#     dplyr::select(slope_df, vegType, bioregion) %>%
#     unnest(slope_df) %>% 
#     dplyr::filter(!is.na(slope)) %>%  # remove the first year (1990)
#     # add slope class column, labels are in percent
#     mutate(
#       slope_class = cut(
#         slope,
#         breaks = c(-Inf, -0.05, -0.02, -0.01, -0.005, 0, 0.005, 0.01, 0.02, 0.05, Inf),
#         labels = c("-5)", "[-5, -2)", "[-2, -1)", "[-1, -0.5)", "[-0.5, 0)",
#                    "[0, 0.5)", "[0.5, 1)", "[1, 2)", "[2, 5)", "[5"),
#         right = FALSE
#       )
#     )%>% 
#     # format for ggplot2
#     group_by(vegType,year,slope_class) %>% 
#     summarise(
#       n = n(),
#       mean_slope = mean(slope, na.rm=T),
#       .groups = "drop") %>% 
#     group_by(vegType,year) %>%
#     mutate(n_total = sum(n),
#            n_ratio = n / n_total*100) %>%
#     ungroup() %>% 
#     mutate(n_ratio=ifelse(mean_slope < 0, -n_ratio, n_ratio)) 
#   return(df_processed)
# }

# Summarise directions of changes at the break points detected by segmented regression
classify_change_direction <- function(tibble_results){
  df_summarised <- tibble_results %>%
    filter(k > 0) %>%
    mutate(slope_transitions = map2(fitted, random_id, function(fitted_data, pid) {
      # get vector of y fitted by segmented regression
      fit_vec <- fitted_data$fit
      names(fit_vec) <- fitted_data$year
      
      # calculate change per year by difference
      slopes <- c(NA, diff(fit_vec))
      
      # choose change points
      change_vec <- abs(diff(slopes)) > 0.001
      bp_positions <- which(change_vec)
      
      # for each bp, categorise change type and return as tibble
      if(length(bp_positions) == 0) return(NULL)
      
      map_dfr(bp_positions, function(bp) {
        change_direction <- ifelse(slopes[bp] < slopes[bp + 1], "inc", "dec")
        before_sign  <- ifelse(slopes[bp] > 0, "pos", "neg")
        after_sign   <- ifelse(slopes[bp + 1] > 0, "pos", "neg")
        change_type <- paste(before_sign, after_sign, change_direction, sep = "_")
        tibble(year = as.numeric(names(slopes[bp])), change_type = change_type)
      })
    })) %>%
    dplyr::select(random_id, vegType, bioregion, slope_transitions) %>% 
    unnest(slope_transitions)
  return(df_summarised)
}

# Summarise and count the types of break points per vegetation type
summarise_change_pattern <- function(df){
  summarised_df <- df %>%
    mutate(change_type_label = recode(change_type,
                                      "pos_pos_dec" = "browining_slowdown",
                                      "pos_pos_inc" = "browining_acceralate",
                                      "pos_neg_dec" = "from_browning_to_greening",
                                      "neg_neg_inc" = "greening_slowdown",
                                      "neg_neg_dec" = "greening_acceralate",
                                      "neg_pos_inc" = "from_greening_to_browning"
    )) %>% 
    group_by(vegType, year,change_type_label) %>% 
    summarise(n = n(),
              .groups = "drop") %>% 
    pivot_wider(names_from=change_type_label,
                values_from=n,
                values_fill=0) %>% 
    pivot_longer(cols=!c(vegType,year),
                 names_to="change_type",
                 values_to="n") %>%
    left_join(df_count, by="vegType") %>% 
    mutate(n_ratio = (n / N)*100)
  return(summarised_df)
}

# Plot fit of segmented regression on sampled pixels
plot_segmented <- function(df, n){
  for(i in 1:n){
    pixel <- df[i,]
    year <- pixel$data[[1]]$plant_year
    y <- pixel$data[[1]]$bareground*100
    y_fit <- pixel$fitted[[1]]$fit*100
    
    plot(year,y, xlab="Year",ylab="Bareground(%)", ylim=c(0,100))
    lines(x=year,y=y_fit, col="red")
    mtext(text=paste0("k=",pixel$k," "),adj=1,line=0.2, cex=0.9)
    mtext(text=pixel$vegType, adj=0,line=1, cex=0.7)
    mtext(text=paste0("S: ",  round(pixel$lat, digits=5),
                      ", E: ",round(pixel$lon, digits=5)),
                      adj=0,line=0.2, cex=0.5)
    
    if(!is.na(unlist(pixel$breakpoints)[[1]])){
      for(yr in pixel$breakpoints){
        abline(v=yr, col="blue", lty="dotted")
        text(x=yr, y=0.05, labels=round(yr), srt=90, adj=c(-0.15,-0.3), cex=0.8)}
    }
  }
}

