######### Functions for pls regression with caret ##############################
# ------------------------------------------------------------------------------
# Function to train a Partial Least Squares (PLS) regression model using caret
# Input:
#   - dataframe: A data frame containing the response and explanatory variables.
#                The response variable must be named "BareGround".
#   - predictors: A character vector of explanatory variable names to use in the model.
#   - blocks: A list of spatial resampling indices, typically from CAST::knndm or similar.
# Output:
#   - A trained caret model object using PLS regression with cross-validation.
func_caret_pls <- function(dataframe, predictors, blocks, n_comp){
  caret_model <- train(
    x = dataframe[, predictors],
    y = dataframe$BareGround,
    method = "pls",                           # Use Partial Least Squares
    preProcess = c('center', 'scale'),        # Standardize predictors before training
    importance = TRUE,                        # Enable variable importance calculation
    metric = "RMSE",
    trControl = trainControl(
      method = "cv",                          # Cross-validation
      # summaryFunction = func_q2,
      selectionFunction = "oneSE",            # Select simplest model within 1-SE of best
      index = blocks$indx_train,               # Cross-validation folds (training)
      indexOut = blocks$indx_test,             # Cross-validation folds (testing)
      savePredictions = "final"               # Save final predictions
    ),
    tuneGrid = data.frame(ncomp = n_comp)
  )
}

# ------------------------------------------------------------------------------
# Function to perform variable selection using PLS and VIP, evaluating model performance
# Input:
#   - list_data: A list containing training and validation data frames, and a spatial blocking object.
#   - predictors: A character vector of predictor variable names.
# Output:
#   - A list containing:
#       1. metrics: A data frame with performance metrics for each VIP-based variable subset.
#       2. vip_full_model: A named numeric vector of VIP scores from the full model.

func_pls_var_selection <- function(list_data, predictors, n_comp, seed=1234) {
  set.seed(seed)
  full_model <- func_caret_pls(list_data$train_df, predictors, list_data$spatial_cv_block, n_comp)
  vip_scores <- VIP(full_model$finalModel, unlist(full_model$bestTune))
  vip_ordered <- sort(vip_scores) %>% names()

  df_metrics <- full_model$results %>%
    mutate(variables = "All")

  for (i in 1:(length(vip_ordered) - 2)) {
    sub_model <- func_caret_pls(
      list_data$train_df,
      vip_ordered[-c(1:i)],
      list_data$spatial_cv_block,
      n_comp
    )

    metrics <- sub_model$results %>%
      mutate(
        variables = paste(sub_model$finalModel$xNames, collapse = ",")
      )

    df_metrics <- bind_rows(df_metrics, metrics)
  }

  return(list(metrics = df_metrics, vip_full_model = vip_scores))
}

# ######### Functions for pls regression with SMOTEr with plsr #################
# ------------------------------------------------------------------------------
# Function to compute standard error
# Input:
#   - x: A numeric vector.
# Output:
#   - Standard error of the input vector.
standard_error <- function(x) {
  sd(x) / sqrt(length(x))
}

# ------------------------------------------------------------------------------
# Function to perform PLS regression with SMOTE-based oversampling and spatial cross-validation
# Input:
#   - dataframe: A data frame containing the training data.
#   - predictors: A character vector of predictor variable names.
#   - data_split: A list containing spatial cross-validation indices.
#   - seed: A numeric value for reproducibility.
# Output:
#   - A list containing:
#       1. results: A data frame summarizing RMSE, MAE, and R2 across components.
#       2. ncomp: Selected number of components based on RMSE + SE rule.
func_pls_smote <- function(dataframe, predictors, data_split, seed) {
  results <- data.frame(fold = integer(),
                        ncomp = numeric(),
                        rmse = numeric(),
                        mae = numeric(),
                        r2 = numeric())
  df_preds <- data.frame("1 comps"=numeric(),
                         "2 comps"=numeric(),
                         "3 comps"=numeric())

  temp_df <- dataframe[, c("BareGround", predictors)]

  for (i in 1:5) {
    train_idx <- data_split$spatial_cv_block$indx_train[[i]]
    test_idx  <- data_split$spatial_cv_block$indx_test[[i]]

    train_data <- temp_df[train_idx, ]
    test_data  <- temp_df[test_idx, ]

    rel_mat <- create_relevance_matrix(train_data$BareGround)

    set.seed(seed)
    smoter_train <- SmoteRegress(BareGround ~ ., dat = train_data,
                                 rel = rel_mat,
                                 C.perc = "extreme",
                                 k = 5)

    pls_model <- plsr(BareGround ~ ., data = smoter_train, scale = TRUE, center = TRUE)

    ncomp_limit <- min(5, length(predictors))
    preds <- predict(pls_model, newdata = test_data)[,,1:ncomp_limit]

    rmse <- apply((preds - test_data$BareGround)^2, 2, mean) %>% sqrt()
    mae  <- apply(abs(preds - test_data$BareGround), 2, mean)
    r2   <- cor(preds, test_data$BareGround)^2

    results <- rbind(results,
                     data.frame(fold = i, ncomp = 1:ncomp_limit,
                                rmse = rmse, mae = mae, r2 = r2))
    df_preds <- rbind(df_preds, preds[,1:3])
  }

  df_res <- results %>%
    group_by(ncomp) %>%
    summarise(
      RMSE = mean(rmse), Rsquared = mean(r2), MAE = mean(mae),
      RMSESD = sd(rmse), RsquaredSD = sd(r2), MAESD = sd(mae),
      .groups = "drop"
    ) %>%
    data.frame()

  min_idx   <- which.min(df_res$RMSE)
  min_rmse  <- df_res$RMSE[min_idx]
  se_rmse   <- standard_error(df_res$RMSE[1:min_idx])
  min_ncomp <- df_res$ncomp[min_idx]

  selected_ncomp <- df_res$ncomp[
    df_res$ncomp <= min_ncomp & df_res$RMSE <= (min_rmse + se_rmse)
  ]

  return(list(results = df_res, ncomp = min(selected_ncomp), preds=df_preds))
}

# ------------------------------------------------------------------------------
# Function to perform variable selection using PLS and VIP, evaluating model performance
# Input:
#   - list_data: A list containing training and validation data frames, and a spatial blocking object.
#   - predictors: A character vector of predictor variable names.
#   - n_comp: The number of components to use in PLS.
#   - seed: A numeric value for reproducibility.
# Output:
#   - A list containing:
#       1. metrics: A data frame with performance metrics for each VIP-based variable subset.
#       2. vip_full_model: A named numeric vector of VIP scores from the full model.
func_pls_smote_var_selection <- function(list_data, predictors, n_comp, seed=1234) {
  temp_df <- list_data$train_df[, c("BareGround", predictors)]
  rel_mat <- create_relevance_matrix(temp_df$BareGround)

  set.seed(seed)
  smoter_train <- SmoteRegress(BareGround ~ ., dat = temp_df,
                               rel = rel_mat,
                               C.perc = "extreme",
                               k = 5)

  full_model <- plsr(BareGround ~ ., data = smoter_train,
                     ncomp = n_comp, scale = TRUE, center = TRUE)

  vip_scores <- VIP(full_model, n_comp)
  vip_ordered <- names(sort(vip_scores))

  full_metrics <- func_pls_smote(temp_df, predictors, list_data, seed = seed)$results[n_comp, ] %>%
    mutate(variables = "All")

  df_metrics <- full_metrics

  for (i in 1:(length(vip_ordered) - 3)) {
    selected_predictors <- vip_ordered[-c(1:i)]

    sub_model <- func_pls_smote(temp_df, selected_predictors, list_data, seed = seed)

    metrics <- sub_model$results[n_comp, ] %>%
      mutate(variables = paste(selected_predictors, collapse = ","))

    df_metrics <- bind_rows(df_metrics, metrics)
  }

  return(list(metrics = df_metrics, vip_full_model = vip_scores))
}
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Function to generate variable selection plot for PLS model
# Inputs:
#   - pls_results: List output from PLS analysis, containing 'metrics' and 'vip_full_model'
#   - label_model: Label for model ID
#   - label_season: Label for season (e.g., Growing/Non-growing)
#   - label_var: Label for predictor type (e.g., "predictors: VI")
#   - label_prep: Label for preprocessing method
# Output:
#   - A ggplot object showing RMSE ± SE vs. number of variables, with 1SE rule selection
# ------------------------------------------------------------------------------
plot_pls_var <- function(pls_results, label_model, label_season, label_var, label_prep) {

  # Extract relevant rows (VIP > 1.0 + intercept)
  df <- pls_results$metrics %>%
    slice(0:(sum(pls_results$vip_full_model <= 1.0) + 1)) %>%
    mutate(n_var = length(vi) - row_number() + 1)

  # 1SE rule threshold
  min_row   <- df[which.min(df$RMSE), ]
  se_min    <- min_row$RMSESD / sqrt(5)  # Assuming 5-fold CV
  threshold <- min_row$RMSE + se_min

  # Select the simplest model within 1SE threshold
  selected <- df %>%
    mutate(selected = RMSE < threshold) %>%
    filter(selected) %>%
    arrange(n_var) %>%
    slice(1) %>%
    dplyr::select(n_var, RMSE)

  # Create plot
  g <- df %>%
    mutate(
      ymin = RMSE - RMSESD / sqrt(5),
      ymax = RMSE + RMSESD / sqrt(5)
    ) %>%
    ggplot(aes(x = n_var, y = RMSE)) +
    geom_linerange(aes(ymin = ymin, ymax = ymax)) +
    geom_point(size = 2) +
    geom_hline(aes(yintercept = threshold), color = "blue", linetype = "dashed") +
    geom_point(data = selected, aes(x = n_var, y = RMSE), color = "red", size = 2) +
    theme_bw() +
    labs(
      y = "RMSE ± SE",
      x = "Number of predictor variables"
    ) +
    ylim(c(14, 21)) +
    annotate("text", x = min(df$n_var), y = Inf, label = label_model, hjust = "left", vjust = 2, size = 4) +
    annotate("text", x = max(df$n_var), y = Inf, label = label_season, hjust = "right", vjust = 2, size = 3) +
    annotate("text", x = max(df$n_var), y = Inf, label = label_var, hjust = "right", vjust = 4, size = 3) +
    annotate("text", x = max(df$n_var), y = Inf, label = label_prep, hjust = "right", vjust = 6, size = 3)

  return(g)
}

