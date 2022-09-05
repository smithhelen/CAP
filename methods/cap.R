#### Prepare test data and training data for CAP method ####

# Load functions
source("methods/helpers.R")       # Functions used internally in the methods

epsilon <- sqrt(.Machine$double.eps)

# Function to map from variable levels to scores.

# Output is both the score information for the variable (i.e. a vector of the same length as the input vector) and
# the level mapping data.frame (mapping from var_level to score)

factor_to_CAP_score <- function(var, dist, class, axes, mp) {
  var_levels <- droplevels(var)
  N <- nlevels(var_levels)
  d.train <- dist[levels(var_levels),levels(var_levels),drop=FALSE]
  A.train <- -0.5 * d.train^2
  B.train <- dbl_center(A.train)  #B is same as Gowers matrix G
  eigen_B <- eigen_decomp(B.train, symmetric=TRUE)
  # use only non-zero eigenvalues
  nlambdas <- sum(eigen_B$values > epsilon)
  if (nlambdas == 0) {
    # No non-zero eigenvectors
    return(NULL)
  }
  ct <- table(Var_level=var_levels,Class=class)
  H <- hat(ct, k=2)
  lambda_B <- filter_eigenvalues(eigen_B$values[seq_len(nlambdas)], axes=axes, mp=mp)
  Q <- eigen_B$vectors[, seq_along(lambda_B), drop=FALSE]
  QHQ <- t(Q) %*% H %*% Q   # Combine Q and H to get Q_score
  eigen_QHQ <- eigen_decomp(QHQ, symmetric=TRUE)
  # also need to drop eigen values here that are too small
  lambda_QHQ <- filter_eigenvalues(eigen_QHQ$values, mp=mp)
  U <- eigen_QHQ$vectors[,seq_along(lambda_QHQ),drop=FALSE]
  Q_score <- Q %*% U
  # Fill Q_scores to individual isolates.  
  score <- left_join(data.frame(Var_level = var_levels), data.frame(Q_score) |> rownames_to_column("Var_level"), by = "Var_level")
  Output <- score |> dplyr::select(-Var_level)
  list(output = Output,
       extra = list(d=dist, var_levels=levels(var_levels), diag.B.train=diag(B.train), lambda_B=lambda_B, 
                    U=U, Q=Q, Q_score=Q_score, num_vars = ncol(Q_score), var_names = colnames(Q_score)))
}

prepare_training_cap <- function(data, vars, class, d, axes, mp) {
  # pull out our var_cols and class
  var_cols <- dplyr::select(data,{{vars}})
  classes   <- data |> pull({{class}})
  # iterate over the var columns and distance matrices, and convert
  prepped <- map2(var_cols, d, factor_to_CAP_score, class = classes, axes, mp) |> compact() # removes empties
  output <- map(prepped,"output")
  prepped_data <- bind_cols(data.frame(classes) |> setNames(data |> select({{class}}) |> colnames()), 
                            map2(output, names(output), ~ .x |> set_names(paste(.y, names(.x), sep="."))))
  extra <- map(prepped, "extra")
  list(training = prepped_data,
       extra = extra)
}

predict_cap <- function(new.var_level, d, diag.B.train, Q, lambda_B, U) {
  # new.var_level is length one
  d.new <- d[new.var_level, names(diag.B.train)]
  d.gower <- diag.B.train - (d.new ^ 2)
  newQ <- d.gower %*% Q / (2 * lambda_B)
  newQscore <- newQ %*% U
  new_var_score <- data.frame(Var_level = new.var_level, newQscore)
  new_var_score
}

# given variable information and a level map, do the mapping. 
impute_score_cap <- function(var, extra) {
  var_levels <- pluck(extra, "var_levels")
  var <- droplevels(var)
  new.var_levels <- setdiff(levels(var), var_levels)
  d <- pluck(extra, "d")
  diag.B.train <- pluck(extra, "diag.B.train")
  lambda_B <- pluck(extra, "lambda_B")
  Q <- pluck(extra, "Q")
  Q_score <- pluck(extra, "Q_score")
  U <- pluck(extra, "U")
  new_scores <- map_df(new.var_levels, ~predict_cap(new.var_level = {.}, d, diag.B.train, Q, lambda_B, U))
  var_level_score <- rbind(data.frame(Q_score) |> rownames_to_column("Var_level"), new_scores)
  list(test_score = data.frame(Var_level = var) |> 
         left_join(var_level_score, by = "Var_level") |>
         select(-Var_level))
}

prepare_test_cap <- function(data, extra, id) {
  id <- data |> select({{id}})
  var_cols <- data |> select(any_of(names(extra)))
  newdata_score <- map2(var_cols, extra, impute_score_cap)
  output <- map(newdata_score,"test_score")
  newdata_pred <- map2_dfc(output, names(output), ~ .x |> set_names(paste(.y, names(.x), sep=".")))
  newdata_pred <- bind_cols(id, newdata_pred)
}

