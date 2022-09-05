#### Prepare test data and training data for pco method ####

# Load functions
source("methods/helpers.R")       # Functions used internally in the methods

epsilon <- sqrt(.Machine$double.eps)

# Function to map from variable levels to scores.

# Output is both the score information for the variable (i.e. a vector of the same length as the input vector) and
# the level mapping data.frame (mapping from var_level to score)
factor_to_pco_score <- function(var, dist, axes) {
  var_levels <- droplevels(var)
  d.train <- dist[levels(var_levels),levels(var_levels),drop=FALSE]
  A.train <- -0.5 * d.train^2
  B.train <- dbl_center(A.train)
  eigen_B <- eigen_decomp(B.train, symmetric=TRUE)
  # use only non-zero eigenvalues
  nlambdas <- sum(eigen_B$values > epsilon)
  if (nlambdas == 0) {
    # No non-zero eigenvectors
    return(NULL)
  }
  # Restrict to a maximum of eigenvectors set by "axes" (default is 2)
  nlambdas <- min(nlambdas, axes)
  # Scale eigenvectors
  lambdas_B <- eigen_B$values[seq_len(nlambdas)]
  X <- sweep(eigen_B$vectors[, seq_len(nlambdas), drop=FALSE], 2, sqrt(abs(lambdas_B)), "*")
  score <- left_join(data.frame(Var_Level = var_levels), data.frame(X) |> rownames_to_column("Var_Level"), by = "Var_Level")
  output <- score |> select(-Var_Level)
  list(output = output,
       extra = list(d=dist, var_levels=levels(var_levels), diag.B.train=diag(B.train), X=X, lambdas_B=lambdas_B, 
                    num_vars = ncol(X), var_names = colnames(X)))
}

prepare_training_pco <- function(data, vars, class, d, axes=2) {
  # pull out our var_cols and class
  var_cols <- dplyr::select(data,{{vars}})
  classes   <- data |> pull({{class}})
  # iterate over the var columns and distance matrices, and convert
  prepped <- map2(var_cols, d, factor_to_pco_score, axes) |> compact() # removes empties
  output <- map(prepped,"output")
  prepped_data <- bind_cols(data.frame(classes) |> setNames(data |> select({{class}}) |> colnames()), 
                            map2(output, names(output), ~ .x |> set_names(paste(.y, names(.x), sep="."))))
  extra <- map(prepped, "extra")
  list(training = prepped_data,
       extra = extra)
}

predict_pco <- function(new.var_level, d, diag.B.train, X, lambdas_B) {
  # new.var_level is length one
  d.new <- d[new.var_level, names(diag.B.train)]
  d.gower <- diag.B.train - (d.new ^ 2)
  newX <- d.gower %*% X / (2 * lambdas_B)
  colnames(newX) = colnames(X)
  new_var_level_score <- data.frame(Var_Level = new.var_level, newX)
  new_var_level_score
}

# given var information and extra information, do the mapping. 
impute_score_pco <- function(var, extra) {
  var <- droplevels(var)
  var_levels <- pluck(extra, "var_levels")
  new.var_levels <- setdiff(levels(var), var_levels)
  d <- pluck(extra, "d")
  diag.B.train <- pluck(extra, "diag.B.train")
  X <- pluck(extra, "X")
  lambdas_B <- pluck(extra, "lambdas_B")
  new_scores <- map_df(new.var_levels, ~predict_pco(new.var_level = {.}, d, diag.B.train, X, lambdas_B))
  var_level_score <- bind_rows(data.frame(X) |> rownames_to_column("Var_Level"), new_scores)
  list(test_score = data.frame(Var_Level = var) |> 
         left_join(var_level_score, by = "Var_Level") |>
         select(-Var_Level))
}

prepare_test_pco <- function(data, extra, id) {
  id <- data |> select({{id}})
  var_cols <- data |> select(any_of(names(extra)))
  newdata_score <- map2(var_cols, extra, impute_score_pco)
  output <- map(newdata_score,"test_score")
  newdata_pred <- map2_dfc(output, names(output), ~ .x |> set_names(paste(.y, names(.x), sep=".")))
  newdata_pred <- bind_cols(id, newdata_pred)
}
