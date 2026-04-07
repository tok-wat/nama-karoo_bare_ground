# ------------------------------------------------------------------------------
# Utility functions

# Making prediction from the saved model with validation dataset
func_pred_valid <- function(list_model, dataframe){
  output <- dataframe %>% 
    dplyr::select(BareGround) %>% 
    rownames_to_column(var="ID") %>%
    mutate(ID=as.character(ID))
  
  for (model in list_model) {
    temp_model <- readRDS(model)
    model_name <- model %>%
      basename() %>%
      str_split("_", simplify = TRUE) %>%
      .[1:2] %>%
      paste(collapse = "_")
    
    if(class(temp_model)=="train"){                         # caret
      temp_df <- dataframe %>% 
        dplyr::select(all_of(c("BareGround",temp_model$finalModel$xNames))) %>% 
        na.omit()
      temp_pred <- predict(temp_model, temp_df)
    } else if (class(temp_model)=="plsRbetamodel") {        # plsR-beta
      temp_df <- dataframe %>% 
        dplyr::select(all_of(c("BareGround",colnames(temp_model$dataX)))) %>% 
        na.omit()
      temp_pred <- predict_plsRbeta(temp_model, temp_df)*100
    } else {
      temp_df <- dataframe %>% 
        dplyr::select(all_of(
          colnames(temp_model$model)
        )) %>% 
        na.omit()
      temp_pred <- predict(temp_model, temp_df)[,,2]  # mvr or pls
    }
    temp_colnames <- c(colnames(output), model_name)
    temp_pred_df <- data.frame(temp = temp_pred) %>% 
      mutate(ID=row.names(temp_df))
    output <- output %>% 
      left_join(temp_pred_df, by="ID")
    colnames(output) <- temp_colnames
  }
  output <- output %>% 
    dplyr::select(-ID)
  return(output)
}

# calculating performance metrics
func_eval_metrics <- function(dataframe){
  dataframe <- na.omit(dataframe)
  df_output <- data.frame(
    model_name=character(),
    r_sq=numeric(),
    rmse=numeric(),
    mae=numeric()
  )
  for(i in 2:ncol(dataframe)){
    bg <- dataframe$BareGround
    temp_lm <-lm(dataframe[,i] ~ bg)
    temp_df <- data.frame(
      model_name=colnames(dataframe)[i],
      r_sq=summary(temp_lm)$r.squared,
      rmse=mean((dataframe[,i]-bg)^2) %>% 
        sqrt() ,
      mae=mean(abs(dataframe[,i]-bg))
    )
    df_output<- rbind(df_output, temp_df)
  }
  return(df_output)
}

# -------------------------------------------------------------------------------
# Visualisation

# Observation vs prediction scatter plot
plot_obs_pred <- function(test_vec, pred_vec, xlab, y_lim=c(-20,120)){
  plot(test_vec, pred_vec,
       xlim =c(0,100),ylim=y_lim,
       xlab = xlab, ylab = "Model Prediction (%)" , col=NULL, bg=rgb(0, 0, 1, alpha=0.3), pch=21, cex=1.5
  )
  temp_lm <-lm(pred_vec ~ test_vec)
  abline(0,1, col="black", lwd=2, lty=2)
  abline(temp_lm, col="red", lwd=2)
  r_sq <- round(summary(temp_lm)$r.squared,3)
  rmse <- mean((test_vec-pred_vec)^2) %>%
    sqrt() %>%
    round(3)
  mae <- mean(abs(test_vec-pred_vec)) %>%
    round(3)
  mtext(
    bquote("RMSE = "~.(rmse)),
    adj=0.05,
    line= -1.8,
    cex = 0.8
  )
  mtext(
    bquote("MAE = "~.(mae)),
    adj=0.05,
    line= -3.0,
    cex = 0.8
  )
  mtext(
    bquote("R"^"2"~"="~.(r_sq)),
    adj=0.05,
    line= -4.2,
    cex = 0.8
  )
}

# Observation vs prediction with landform category for transect data
plot_slope_flat <- function(dataframe, valid_vec, xlab, y_lim = c(-20,120)){
  df_temp <- dataframe %>%
    dplyr::select(BareGround,Notes) %>%
    mutate(pred = valid_vec) %>%
    na.omit()
  df_flat_temp <- df_temp %>%
    filter(str_detect(Notes, "Plain"))
  df_slope_temp <- df_temp %>%
    filter(str_detect(Notes, "Slope"))
  df_others_temp <- df_temp %>%
    filter(!str_detect(Notes, "Slope")) %>%
    filter(!str_detect(Notes, "Plain"))

  # plot flat points with blue circles
  plot(df_flat_temp$BareGround, df_flat_temp$pred,
       xlim =c(0,100),ylim=y_lim,
       xlab = xlab, ylab = "Model Prediction (%)" , col=NULL, bg=rgb(0, 0.8, 0.4, alpha=0.6), pch=21, cex=1.5
  )
  # plot slope points with red triangles
  points(df_slope_temp$BareGround, df_slope_temp$pred, col=NULL, bg=rgb(1, 0, 0, alpha=0.6), pch=24, cex=1.5)
  points(df_others_temp$BareGround, df_others_temp$pred, col=NULL, bg=rgb(0, 0, 0, alpha=0.6), pch=22, cex=1.5)

  temp_lm <-lm(df_temp$pred ~ df_temp$BareGround)
  abline(0,1, col="black", lwd=2, lty=2)
  abline(temp_lm, col="red", lwd=2)
  r_sq <- round(summary(temp_lm)$r.squared,3)
  rmse <- mean((df_temp$pred-df_temp$BareGround)^2) %>%
    sqrt() %>%
    round(3)
  mae <- mean(abs(df_temp$pred-df_temp$BareGround)) %>%
    round(3)
  mtext(
    bquote("RMSE = "~.(rmse)),
    adj=0.05,
    line= -1.8,
    cex = 0.8
  )
  mtext(
    bquote("MAE = "~.(mae)),
    adj=0.05,
    line= -3.0,
    cex = 0.8
  )
  mtext(
    bquote("R"^"2"~"="~.(r_sq)),
    adj=0.05,
    line= -4.2,
    cex = 0.8
  )
}

# -----------------------------------------------------------------------
# ANOVA to test the effects of preprocess, season and algorithms on
# the performance metrics
# Table S3

# ANOVA with lm framework
func_anova_summary <- function(lm_model){
  lm_est <- cbind(coef(summary(lm_model))[-1,1:2]) %>%
    as.data.frame() %>% 
    mutate(dep_var =  c("smote", "box_cox", "season","Algo_plsbeta", "Algo_rf")) 
  
  Anova(lm_model, type=2) %>% 
    data.frame() %>% 
    rownames_to_column(var="dep_var") %>% 
    mutate(eta_sq = Sum.Sq/sum(Sum.Sq)) %>% 
    full_join(lm_est, by="dep_var") %>% 
    slice(c(1:4,6:7,5)) %>%  # recorder the column
    dplyr::select(dep_var, Df,Estimate, `Std. Error`, F.value, Pr..F.,eta_sq) %>% 
    return()
}

func_model_comparison <- function(df){
  # Interaction not considered as the models were insignificant from 
  # the models without interaction terms
  lm_rmse <- lm(rmse ~ smote + box_cox + season + Algorithm, data=df)
  lm_mae  <- lm(mae  ~ smote + box_cox + season + Algorithm, data=df)
  lm_r_sq <- lm(r_sq ~ smote + box_cox + season + Algorithm, data=df)
  
  out_df <- cbind(func_anova_summary(lm_rmse), 
                  func_anova_summary(lm_mae)[,-(1:2)],
                  func_anova_summary(lm_r_sq)[,-(1:2)])
  
  colnames(out_df) <- c("dep_var", "df", 
                        "RMSE_estimate", "RMSE_se", "RMSE_F_value", "RMSE_p_value", "RMSE_eta_sq",
                        "MAE_estimate",  "MAE_se",  "MAE_F_value",  "MAE_p_value",  "MAE_eta_sq",
                        "R2_estimate",   "R2_se",   "R2_F_value",   "R2_p_value",   "R2_eta_sq")
  return(out_df)
}

# Tukey pairwise test among algorithms
func_posthoc_tukey <- function(df){
  tk_rmse <- aov(scale(rmse) ~ smote + box_cox + season + Algorithm, data=df) %>% 
    TukeyHSD()
  tk_mae  <- aov(scale(mae) ~ smote + box_cox + season + Algorithm, data=df) %>% 
    TukeyHSD()
  tk_r_sq <- aov(scale(r_sq) ~ smote + box_cox + season + Algorithm, data=df) %>% 
    TukeyHSD()
  rbind(tk_rmse$Algorithm,tk_mae$Algorithm,tk_r_sq$Algorithm) %>% 
    as.data.frame() %>% 
    mutate(metrics = rep(c("RMSE","MAE","R²"), times=1, each=3),
           pairs   = rep(c("PLS - PLS-beta","PLS - RF","PLS-beta - RF"), times=3, each=1)) %>%
    arrange(desc(pairs)) %>% 
    return()
} 

# visualisation of Tukey's test
gg_posthoc <- function(df){
  df %>% 
    ggplot(aes(x = x_pos, y = diff)) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width=0.2) +
    geom_point(size = 3) +
    geom_hline(aes(yintercept = 0), color = "red") +
    coord_flip() +     
    scale_x_continuous(
      breaks = df$x_pos,
      labels = df$label
    )+
    geom_text(aes(x_pos, 2.7,
                  label = paste0(sprintf("%.2f", diff),
                                 " (",
                                 sprintf("%.2f", lwr),
                                 ", ",
                                 sprintf("%.2f", upr),
                                 ")")
    ),size =3.3
    )+
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank()
    )+
    labs(y = "Standardised difference in mean levels and 95% CI",
         x = "") +
    ylim(c(-3,3))
}

# -----------------------------------------------------------------------
# ANCOVA to test the effect of landforms
# Table S4

func_ancova <- function(df){
  output <- list()
  
  for(i in 3:(ncol(df)-1)){
    temp_df <- df[,c(2,i,ncol(df))] %>% 
      na.omit()
    colnames(temp_df) <- c("BareGround", "temp_x", "Landforms")
    
    #temp_lm <- lm(BareGround ~ temp_x + Landforms, data=temp_df)
    temp_lm <- lm(BareGround ~ temp_x * Landforms, data=temp_df)
    #p_interaction <- anova(temp_lm_interaction)$`Pr(>F)`[3]
    anova_type2 <- Anova(temp_lm, type = 2) %>% 
      data.frame() %>% 
      rownames_to_column() %>% 
      mutate(model = colnames(df[i]),
             beta = c(coef(summary(temp_lm))[-1,1], NA),
             se   = c(coef(summary(temp_lm))[-1,2], NA),
             # p_int = p_interaction,
             partial_eta2 = c(Sum.Sq[1]/(Sum.Sq[1]+Sum.Sq[4]),  # =r2 of lm model
                              Sum.Sq[2]/(Sum.Sq[2]+Sum.Sq[4]),
                              Sum.Sq[3]/(Sum.Sq[3]+Sum.Sq[4]),
                              NA)) 
    output[[i-2]] <- anova_type2
  }
  
  # formatting
  output <- list_rbind(output) %>% 
    dplyr::select(model, rowname, beta, se, Df, F.value, Pr..F., Sum.Sq, partial_eta2) %>% 
    rename(dep_var = "rowname") %>%
    mutate(dep_var = case_when(
      dep_var == "temp_x"   ~ "Predicted value",
      dep_var == "Landforms"   ~ "Landforms (slope)",
      dep_var == "temp_x:Landforms"   ~ "Predicted value:Landforms",
      .default = as.character(dep_var)
    ))
  
  return(output)
}
