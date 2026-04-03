######### Functions for randomforest_regression ################################
# ------------------------------------------------------------------------------
# Function to select predictors based on Variance Inflation Factor (VIF)
# Input:
#   - dataframe: A data frame containing the response and explanatory variables.
#                The response variable must be named "BareGround".
#   - predictors: A character vector specifying the names of candidate explanatory variables.
# Output:
#   - A character vector of selected predictors with VIF less than the specified threshold.
func_colinear <- function(dataframe, predictors) {
  pref_order <- dataframe %>%
    dplyr::select(BareGround,all_of(vi)) %>%
    cor() %>%
    as.data.frame() %>%
    dplyr::select(BareGround) %>%
    arrange(desc(abs(BareGround))) %>%
    rownames()
  cor_select(
    df = dataframe,
    predictors = predictors,
    preference_order = pref_order[2:6], # prioritise highly correlated variables
    max_cor = 0.90,
    quiet = TRUE
  ) %>%
    as.vector()
}

# ------------------------------------------------------------------------------
# Function to perform Recursive Feature Elimination (RFE) using Random Forest
# Input:
#   - dataframe: A data frame containing the response and explanatory variables.
#                The response variable must be named "BareGround".
#   - predictors: A character vector specifying the names of candidate explanatory variables.
#   - blocks: A list of spatial resampling indices created with CAST::knndm (or similar).
#   - seed: An integer to set the random seed for reproducibility.
# Output:
#   - A list containing:
#     [[1]]: A data frame of model performance metrics (e.g., RMSE, R²) for each subset size.
#     [[2]]: The number of predictors in the most parsimonious model within one standard error (1-SE) of the best model.
#     [[3]]: A character vector of selected predictors from the 1-SE model.
func_rfe <- function(dataframe, predictors, blocks, seed = 1234) {
  set.seed(seed)  # Ensure reproducibility

  output <- rfe(
    x = dataframe[, predictors],
    y = dataframe$BareGround,
    sizes = 1:length(predictors),
    rfeControl = rfeControl(
      functions = rfFuncs,
      method = "cv",
      index = blocks$indx_train,
      indexOut = blocks$indx_test,
    )
  )

  nparam <- oneSE(output$results, metric = "RMSE", maximize = FALSE, num=5)
  selected_pred <- output$optVariables[1:nparam]

  return(list(
    results = output$results,
    n_selected = nparam,
    selected_predictors = selected_pred
  ))
}

# ------------------------------------------------------------------------------
# Function to train a Random Forest model using the caret package
# Input:
#   - dataframe: A data frame containing the response and explanatory variables.
#                The response variable must be named "BareGround".
#   - predictors: A character vector specifying the names of explanatory variables to use.
#   - blocks: A list of spatial resampling indices created with CAST::knndm (or similar).
#   - seed: An integer to set the random seed for reproducibility.
# Output:
#   - A caret model object trained using cross-validation.
func_caret_rf <- function(dataframe, predictors, blocks, seed = 1234) {
  set.seed(seed)  # Ensure reproducibility

  caret_model <- train(
    x = dataframe[, predictors],
    y = dataframe$BareGround,
    method = "rf",
    # preProcess = c("center", "scale"),
    importance = TRUE,
    trControl = trainControl(
      method = "cv",
      index = blocks$indx_train,
      indexOut = blocks$indx_test,
      savePredictions = "final"
    )
  )

  return(caret_model)
}

# ------------------------------------------------------------------------------
# Function to generate RFE variable selection plot
# Inputs:
#   - rfe_obj: Object returned by RFE, containing $results with RMSE and RMSESD
#   - label_model: Label for model ID
#   - label_season: Label for season (e.g., Growing/Non-growing)
#   - label_var: Label for predictor type
#   - label_prep: Label for preprocessing method
# Output:
#   - A ggplot object showing RMSE ± SE vs. number of variables
# ------------------------------------------------------------------------------
plot_rfe <- function(rfe_obj, label_model, label_season, label_var, label_prep) {
  df <- rfe_obj$results
  
  # 1SE rule threshold
  min_row   <- df[which.min(df$RMSE), ]
  se_min    <- min_row$RMSESD / sqrt(5)  # Assuming 5-fold CV
  threshold <- min_row$RMSE + se_min
  
  selected <- df %>%
    mutate(selected = RMSE < threshold) %>%
    filter(selected) %>%
    arrange(Variables) %>%
    slice(1) %>%
    dplyr::select(Variables, RMSE)
  
  # Create plot
  df %>%
    mutate(
      ymin = RMSE - RMSESD / sqrt(5),
      ymax = RMSE + RMSESD / sqrt(5)
    ) %>%
    ggplot(aes(x = Variables, y = RMSE)) +
    geom_linerange(aes(ymin = ymin, ymax = ymax)) +
    geom_point(size = 2) +
    geom_point(data = selected, aes(x = Variables, y = RMSE), color = "red", size = 2) +
    geom_hline(aes(yintercept = threshold), color = "blue", linetype = "dashed") +
    theme_bw() +
    labs(y = "RMSE ± SE", x = "Number of predictor variables") +
    ylim(c(10, 25)) +
    annotate("text", x = max(df$Variables), y = Inf, label = label_model, hjust = "right", vjust = 2, size = 4) +
    annotate("text", x = max(df$Variables), y = Inf, label = label_season, hjust = "right", vjust = 4.5, size = 3) +
    annotate("text", x = max(df$Variables), y = Inf, label = label_var, hjust = "right", vjust = 6.5, size = 3) +
    annotate("text", x = max(df$Variables), y = Inf, label = label_prep, hjust = "right", vjust = 8.5, size = 3)
}

# ------------------------------------------------------------------------------
# Function to run the Random Forest model pipelinet
# Inputs:
#   - dataframe: a data frame for training dataset
#   - predictors: predictor variables as a vector
#   - blocks: spatial blocking, output of UBL::knndm()
#   - seed: random seeds
# Output:
#   - A list of [[1]] selected variables with a threshold of r <0.9, 
#               [[2]] selected variables with recursive feature elimination,
#               [[3]] developed RF model
# ------------------------------------------------------------------------------
func_rf_pipeline <- function(dataframe, predictors, blocks, seed){
  selected_vars <- func_colinear(dataframe, predictors)
  temp_rfe <- func_rfe(
    dataframe = dataframe,
    predictors = selected_vars, 
    blocks = blocks, 
    seed=seed)
  rf_model <- func_caret_rf(
    dataframe = dataframe,
    predictors = temp_rfe$selected_predictors, 
    blocks = blocks, 
    seed=seed)
  return(list(cor_output = selected_vars,
              RFE_output = temp_rfe, 
              best_model = rf_model))
}

# ------------------------------------------------------------------------------
# Function to perform the full random forest pipeline with SMOTE and feature selection
# Inputs:
#   - dataframe: a data frame for training dataset
#   - predictors: predictor variables as a vector
#   - blocks: spatial blocking, output of UBL::knndm()
#   - seed: random seeds
# Output:
#   - A list of [[1]] selected variables with a threshold of r <0.9, 
#               [[2]] selected variables with recursive feature elimination,
#               [[3]] cross-validation output
#               [[4]] developed RF model
# ------------------------------------------------------------------------------
func_rf_smote_pipeline <- function(dataframe, predictors, blocks, seed = 1234) {
  
  # Step 1: Remove collinear variables
  selected_vars <- func_colinear(dataframe, predictors)
  
  # Step 2: Apply SMOTE with spatial cross-validation blocks
  temp_smote <- func_smote_df(dataframe, selected_vars, blocks, seed)
  
  # Step 3: Perform Recursive Feature Elimination (RFE) for feature selection
  temp_rfe <- func_rfe(
    dataframe = temp_smote$train_df,
    predictors = selected_vars,
    blocks = blocks$spatial_cv_block,
    seed = seed
  )
  
  # Step 4: Prepare the dataset for modeling with selected predictors
  temp_df <- temp_smote$train_df[, c("BareGround", temp_rfe$selected_predictors)]
  
  # Step 5: Cross-validation evaluation
  results <- data.frame(fold = integer(),
                        rmse = numeric(),
                        mae = numeric(),
                        r2 = numeric())
  
  for (i in 1:5) {
    set.seed(seed)
    train_idx <- temp_smote$spatial_cv_block$indx_train[[i]]
    test_idx  <- temp_smote$spatial_cv_block$indx_test[[i]]
    
    train_data <- temp_df[train_idx, ]
    test_data  <- temp_df[test_idx, ]
    
    # Create relevance matrix for SMOTE (even if not used here, may be used in upstream functions)
    rel_mat <- create_relevance_matrix(train_data$BareGround)
    
    # Train Random Forest model
    cv_rf_model <- train(
      x = train_data[, temp_rfe$selected_predictors],
      y = train_data$BareGround,
      method = "rf",
      # preProcess = c("center", "scale"),
      importance = TRUE
    )
    
    # Make predictions
    preds <- predict(cv_rf_model, newdata = test_data)
    
    # Calculate evaluation metrics
    rmse <- sqrt(mean((preds - test_data$BareGround)^2))
    mae  <- mean(abs(preds - test_data$BareGround))
    r2   <- cor(preds, test_data$BareGround)^2
    
    # Store results
    results <- rbind(results,
                     data.frame(fold = i, rmse = rmse, mae = mae, r2 = r2))
  }
  
  # Step 6: Summarize cross-validation results
  df_res <- results %>%
    summarise(
      RMSE = mean(rmse), Rsquared = mean(r2), MAE = mean(mae),
      RMSESD = sd(rmse), RsquaredSD = sd(r2), MAESD = sd(mae),
      .groups = "drop"
    ) %>%
    data.frame()
  
  # Step 7: Train final model on all data with SMOTE
  all_rel_mat <- create_relevance_matrix(dataframe$BareGround)
  
  set.seed(seed)
  all_smote <- SmoteRegress(
    BareGround ~ .,
    dat = dataframe[, c("BareGround", temp_rfe$selected_predictors)],
    rel = all_rel_mat,
    C.perc = "extreme",
    k = 5
  )
  
  set.seed(seed)
  all_rf_model <- train(
    x = all_smote[, temp_rfe$selected_predictors],
    y = all_smote$BareGround,
    method = "rf",
    # preProcess = c("center", "scale"),
    importance = TRUE
  )
  
  # Return all outputs
  return(list(
    cor_output = selected_vars,
    RFE_output = temp_rfe,
    cv_output = df_res,
    best_model = all_rf_model
  ))
}