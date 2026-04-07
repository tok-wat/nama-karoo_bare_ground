######### Variable Names  ######################################################
# Define predictor variable groups
vi <- c("bsi","evi","exg","msavi","ndmi","ndvi","osavi","satvi","savi","swirRatio","tc_bright","tc_green","tc_wet")
minerals <- c("calcite","chlorite","dolomite","goethite","hematite","illite_muscovite","kaolinite") # excluding gypsum and vermiculite
soil <- c("soil_depth","clay","sand","silt","stone", "nitrogen","phosphorous","potassium","ph","soc")
topography <- c("elevation", "slope", "aspect") 


######### Functions for preprocessing  #########################################
# ------------------------------------------------------------------------------
# Function to get parameters for box-cox transformation as a preprocess
# Input:
#   - x: A vector of vegetation index to be used as predictor variables 
#       to be applied box-cox transformation.
#   - y: A vector of 'BareGround' to be used as a response variable.
# Output:
#   - A list with two elements:
#       1. lambda: lambda for box-cox transformation.
#       2. shift_x: value to be added to account for negative value in x.
get_boxcox_params <- function(x, y) {
  shift_x <- ifelse(min(x, na.rm = TRUE) <= 0, abs(min(x, na.rm = TRUE)) + 1e-4, 0)
  x_shifted <- x + shift_x
  temp_y <- y + 1e-4
  bx <- boxcox(temp_y ~ x_shifted, lambda = seq(-0.5, 1.5, len = 100))
  lambda <- bx$x[which.max(bx$y)]
  list(lambda = lambda, shift_x = shift_x)
}

# ------------------------------------------------------------------------------
# Function apply box-cox transformation using parameters estimated
#  with 'get_boxcox_params'
# Input:
#   - x: A vector of vegetation index to be used as predictor variables
#       to be applied box-cox transformation.
#   - lambda: lambda for box-cox transformation. Estimated in 'get_boxcox_params.
#   - shift_x: value to be added to account for negative value in x.
#       Estimated in 'get_boxcox_params.
# Output:
#   - A vector with transformed values
apply_boxcox_transform <- function(x, lambda, shift_x) {
  x_shifted <- x + shift_x
  if (abs(lambda) < 1e-6) {
    transformed <- log(x_shifted)
  } else {
    transformed <- (x_shifted ^ lambda - 1) / lambda
  }
  return(transformed)
}

# ------------------------------------------------------------------------------
# Function to split spatial data into training and validation sets with stratified
# sampling (5-bins) and then create spatial blocks
# Input:
#   - sf_object: An sf object containing spatial point data.
#   - ratio: Proportion of data to assign to training set (default is 0.7).
#   - k: Number of nearest neighbors used in spatial blocking (default is 5).
#   - seed: Random seed for reproducibility (default is 1234).
# Output:
#   - A list containing:
#       1. train_df: A data frame of the training data (geometry dropped).
#       2. valid_df: A data frame of the validation data (geometry dropped).
#       3. spatial_cv_block: A spatial blocking object for cross-validation.
func_data_split <- function(sf_object, ratio = 0.7, k = 5, seed = 1234) {
  set.seed(seed)
  
  sf_bin <- sf_object %>%
    mutate(y_bin = cut(sf_object$BareGround,
                       breaks =quantile(sf_object$BareGround, probs = seq(0, 1, 0.2), na.rm = TRUE),
                       include.lowest=TRUE,
                       labels=FALSE))
  train_idx <- createDataPartition(sf_bin$y_bin, p=ratio, list=FALSE)
  
  sf_train <- sf_object[train_idx,]
  df_train <- st_drop_geometry(sf_train)
  df_valid <- sf_object[-train_idx,] %>%
    st_drop_geometry()
  
  spatial_blocks <- knndm(sf_train, study_area, k = k)
  
  return(list(train_df = df_train, valid_df = df_valid, spatial_cv_block = spatial_blocks))
}

# ------------------------------------------------------------------------------
# Function to create a list of dataframe with oversampling data with SMOTE and
# and spatial blocking
# Input:
#   - dataframe: A data frame containing the response and explanatory variables.
#                The response variable must be named "BareGround".
#   - predictors:  A character vector of predictor variable names.
#   - data_split: A list containing spatial cross-validation indices.
#   - seed: A numeric value for reproducibility.
# Output:
#   - A list containing:
#       1.train_df: A dataframe with oversampled by SMOTE
#       2.spatial_cv_block: A spatial block for spatial cross-validation
func_smote_df <- function(dataframe, predictors, data_split, seed=1234){
  temp_df <- dataframe[,c("BareGround",predictors)]
  
  train_df <- temp_df[0,]
  test_df <- temp_df[0,]
  
  for(i in 1:5){
    train_idx <- data_split$spatial_cv_block$indx_train[[i]]
    test_idx <- data_split$spatial_cv_block$indx_test[[i]]
    
    train_data <- temp_df[train_idx, ]
    test_data  <- temp_df[test_idx, ]
    
    rel_mat <- create_relevance_matrix(train_data$BareGround)
    set.seed(seed)
    smoter_train <- SmoteRegress(BareGround ~ ., dat = train_data,
                                 rel = rel_mat,
                                 C.perc = "extreme",
                                 k = 5)
    temp_train_idx <- 1:nrow(smoter_train)
    temp_train_df <- smoter_train %>%
      mutate(idx= temp_train_idx) %>%
      mutate(block=i) %>%
      mutate(type="train")
    
    temp_test_idx <- 1:nrow(test_data)
    temp_test_df <- test_data %>%
      mutate(idx= temp_test_idx) %>%
      mutate(block=i)%>%
      mutate(type="test")
    
    train_df <- rbind(train_df, temp_train_df)
    test_df <- rbind(test_df, temp_test_df)
  }
  
  smote_df <- rbind(train_df, test_df)%>%
    mutate(idx = row_number())
  train_list <- list()
  test_list <- list()
  
  for(j in 1:5){
    temp_train_id <- smote_df %>%
      filter(type=="train") %>%
      filter(block==j) %>%
      dplyr::select(idx) %>%
      unlist() %>%
      unname()
    temp_test_id <- smote_df %>%
      filter(type=="test") %>%
      filter(block==j) %>%
      dplyr::select(idx) %>%
      unlist() %>%
      unname()
    
    train_list[[j]]<- temp_train_id
    test_list[[j]]<- temp_test_id
  }
  blocks <- list(indx_train=train_list, indx_test=test_list)
  
  smote_df <- smote_df %>%
    dplyr::select(-c(idx,block,type))
  return(list(train_df=smote_df, spatial_cv_block=blocks))
}

# ------------------------------------------------------------------------------
# Function to create a relevance matrix for SMOTER
# Input:
#   - response: A numeric vector of the response variable (e.g., BareGround).
# Output:
#   - A 3x3 matrix defining the relevance of extreme and mean values.
create_relevance_matrix <- function(response) {
  matrix(c(
    0, 1, 0,                         # Minimum extreme value (relevance = 1)
    mean(response), 0, 0,           # Mean value (relevance = 0)
    max(response), 1, 0                       # Maximum extreme value (relevance = 1)
  ), ncol = 3, byrow = TRUE)
}


# ------------------------------------------------------------------------------
# Function to compute predictions from a fitted PLS-beta model
# Input:
#   - plsRbeta_model: A fitted PLS-beta model object
#   - dataframe: A data frame containing predictors for prediction
# Output:
#   - A numeric vector of predicted values (on response scale)
predict_plsRbeta <- function(plsRbeta_model, dataframe){
  temp_df <- dataframe %>%
    mutate(intercept = 1) %>%
    dplyr::select(intercept, all_of(rownames(plsRbeta_model$Coeffs)[-1]))
  
  # Compute predictions using the PLS-beta model
  logit_values <- as.matrix(temp_df) %*% as.vector(plsRbeta_model$Coeffs)
  predictions <- exp(logit_values) / (1 + exp(logit_values))
  return(predictions)
}
