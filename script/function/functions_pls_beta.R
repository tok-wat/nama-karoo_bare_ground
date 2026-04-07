# ######### Functions for pls_beta_regression ##################################
# ------------------------------------------------------------------------------
# Function to select predictor variables with consistent sign in confidence intervals
# Input:
#   - ci: A data frame or matrix containing confidence interval values,
#         where columns 7 and 8 represent the lower and upper bounds.
# Output:
#   - A character vector of selected predictor variable names
get_selected_var <- function(ci){
  # Identify rows where both lower and upper bounds are either positive or negative
  consistent_sign <- (ci[, 7] < 0 & ci[, 8] < 0) | (ci[, 7] > 0 & ci[, 8] > 0)
  # Extract variable names where the condition is TRUE
  selected_vars <- names(consistent_sign[consistent_sign])
  return(selected_vars)
}

# ------------------------------------------------------------------------------
# Function to perform PLS-beta regression with cross-validation folds
# Input:
#   - y: Response variable (numeric vector)
#   - X: Predictor matrix or data frame
#   - fold_list: A list containing train/test indices for CV
#   - ncomp: Number of PLS components
# Output:
#   - A named vector with average RMSE, R², MAE and their standard deviations
func_plsRbeta_fold <- function(y, X, fold_list, ncomp) {
  require(Metrics)

  # output vecotrs
  rmse_list <- numeric(length(fold_list$indx_train))
  mae_list <- numeric(length(fold_list$indx_train))
  r2_list <- numeric(length(fold_list$indx_train))

  for (i in seq_along(fold_list$indx_train)) {
    test_idx <- fold_list$indx_test[[i]]
    train_idx <- fold_list$indx_train[[i]]

    X_train <- X[train_idx, ]
    y_train <- y[train_idx]
    X_test <- X[test_idx, ]
    y_test <- y[test_idx]

    df_test <- X_test %>%
      as.data.frame() %>%
      mutate(BareGround = y_test)

    # model construction and prediction
    model <- PLS_beta(
      dataY=y_train,
      dataX=X_train,
      nt = ncomp,
      modele = "pls-beta",
      link.phi= "log",
      scaleX = TRUE,
      verbose = FALSE)

    y_pred <- predict_plsRbeta(model, df_test)*100
    y_test <- y_test*100

    # calculate metrics
    rmse_list[i] <- rmse(y_test, y_pred)
    mae_list[i] <- mae(y_test, y_pred)
    r2_list[i] <- cor(y_test, y_pred)^2 # to align with R2 calculated in 'caret'

  }

    return(c(
    RMSE = mean(rmse_list),
    Rsquared = mean(r2_list),
    MAE = mean(mae_list),
    RMSESD = sd(rmse_list),
    RsquaredSD = sd(r2_list),
    MAESD = sd(mae_list)
  ))
}

# ------------------------------------------------------------------------------
# Function to evaluate PLS-beta performance over multiple components via CV
# Input:
#   - y: Response vector
#   - X: Predictor matrix/data frame
#   - fold_list: Fold indices for cross-validation
#   - n_comp: Maximum number of components to test
# Output:
#   - A data frame with metrics for each number of components
func_plsRbeta_cv <- function(y, X, fold_list, n_comp=5){
  results <- data.frame(ncomp = numeric(),
                        RMSE = numeric(),
                        Rsquared = numeric(),
                        MAE = numeric(),
                        RMSESD = numeric(),
                        RsquaredSD = numeric(),
                        MAESD = numeric())

  for(i in 1:n_comp){
    cv_model<- func_plsRbeta_fold(y, X, fold_list, i)
    cv_model <- cv_model %>%
      t() %>%
      as.data.frame() %>%
      mutate(ncomp=i)
    results <- rbind(results,cv_model)
  }
  return(results)
}

# ------------------------------------------------------------------------------
# Function to evaluate full model fit with AIC/BIC after CV for component selection
# Input:
#   - dataframe: A data frame with response and predictors
#   - predictors: Character vector of predictor column names
#   - data_split: A list including cross-validation folds
#   - seed: Random seed (not used here but included for compatibility)
# Output:
#   - A data frame with CV results plus AIC and BIC
func_plsRbeta_ncomp_cv<- function(dataframe, predictors, data_split, seed=1234){
  temp_cv <- func_plsRbeta_cv(
    y= dataframe$BareGround,
    X= dataframe[,predictors],
    fold_list=data_split$spatial_cv_block,
    n_comp=5)
  full_model <- PLS_beta(
    dataY= dataframe$BareGround,
    dataX= dataframe[,predictors],
    modele = "pls-beta",
    link.phi= "log",
    nt=5,
    scaleX = TRUE,
    verbose = FALSE
  )

  temp_cv$AIC <- as.vector(full_model$AIC)[-1]
  temp_cv$BIC <- as.vector(full_model$BIC)[-1]

  return(temp_cv)
}

# ------------------------------------------------------------------------------
# Function to evaluate PLS-beta performance using SMOTE-augmented training data
# Input:
#   - dataframe: Original data frame
#   - predictors: Predictor variable names
#   - data_split: List with fold indices for CV
#   - seed: Random seed for SMOTE
# Output:
#   - A data frame with CV results and AIC/BIC for SMOTE-augmented model
func_plsRbeta_smote_ncomp_cv <- function(dataframe, predictors, data_split, seed=1234){
  temp_df_smote <- func_smote_df(
    dataframe=dataframe,
    predictors=predictors,
    data_split=data_split,
    seed=1234)
  temp_smote_cv <- func_plsRbeta_cv(
    y= temp_df_smote$train_df$BareGround,
    X= temp_df_smote$train_df[,predictors],
    fold_list=temp_df_smote$spatial_cv_block,
    n_comp=5)
  smoter_all <- SmoteRegress(
    BareGround ~ ., dat = dataframe[,c("BareGround",predictors)],
    rel = create_relevance_matrix(dataframe$BareGround),
    C.perc = "extreme",
    k = 5)
  full_model <- PLS_beta(
    dataY= smoter_all$BareGround,
    dataX= smoter_all[,predictors],
    modele = "pls-beta",
    link.phi= "log",
    nt=5,
    scaleX = TRUE,
    verbose = FALSE
  )

  temp_smote_cv$AIC <- as.vector(full_model$AIC)[-1]
  temp_smote_cv$BIC <- as.vector(full_model$BIC)[-1]

  return(temp_smote_cv)
}
# ------------------------------------------------------------------------------
# Function to plot cross-validation metrics (AIC, BIC, RMSE, R², MAE)
# across different numbers of PLS components
# Input:
#   - dataframe: A data frame containing metrics for 1 to 5 components
# Output:
#   - Five of line plots
plot_cv <- function(dataframe){
  plot(1:5, dataframe$AIC, type="o",lwd=4, xlab="Number of Components", ylab="AIC",col="red")
  plot(1:5, dataframe$BIC, type="o",lwd=4, xlab="Number of Components", ylab="BIC",col="red")
  plot(1:5, dataframe$RMSE, type="o",lwd=4, xlab="Number of Components", ylab="RMSE",col="blue")
  plot(1:5, dataframe$Rsquared, type="o",lwd=4, xlab="Number of Components", ylab=bquote("R"^"2"),col="blue")
  plot(1:5, dataframe$MAE, type="o",lwd=4, xlab="Number of Components", ylab="MAE",col="blue")
}
