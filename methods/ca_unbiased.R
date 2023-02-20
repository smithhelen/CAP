#### Prepare test data and training data for ca_unbiased method ####

# Function to map from variable levels to ranks.
epsilon <- sqrt(.Machine$double.eps)

# Output is both the score information for the variable (i.e. a vector of the same length as the input vector) and
# the level mapping data.frame (mapping from var_level to score)
factor_to_ca0_score <- function(var, class, axes) {
  var_levels <- droplevels(var)
  if(nlevels(var_levels) < 2){return(NULL)}
  ct <- table(Var_Level=var_levels,Class=class)
  # add zero row - can't add zeros to ct as P will complain. So add 1/num.classes to P (equal probabilities across classes)
  new <- matrix(rep((1/nlevels(class)),times=nlevels(class)),nrow=1, 
                ncol=nlevels(class), dimnames = list(Var_Level="new",Class=colnames(ct)))
  P <- rbind(ct/rowSums(ct),new)
  S <- cov.wt(P, wt = c(rowSums(ct),0))$cov
  eigen_S <- eigen_decomp(S, symmetric=TRUE) ## PCA of weighted covariance matrix of class probabilites
  # Restrict to a maximum of eigenvectors set by "axes" (default is 2)
  nlambdas <- min(sum(eigen_S$values > epsilon), axes)
  # principal components
  pc <- eigen_S$vectors
  X <- P %*% pc[, seq_len(nlambdas), drop=FALSE] |> as.data.frame() 
  score <- left_join(data.frame(Var_Level = var_levels), X |> rownames_to_column("Var_Level"), by = "Var_Level")
  output <- score |> select(-Var_Level)
  list(output = output,
       extra = list(X=X, var_levels=levels(var_levels), num_vars = ncol(X), var_names = colnames(X)))
}

# iterate over the variable columns and grab the output as a new data.frame to send into ca, and store the absent level stuff for later
prepare_training_ca0 <- function(data, vars, class, axes=2, residualised=NULL) {
  # pull out our var_cols and class
  var_cols <- data |> select(all_of(vars))
  classes   <- data |> pull(class)
  # iterate over the var columns, and convert
  prepped <- map(var_cols, factor_to_ca0_score, class = classes, axes) |> compact()
  output <- map(prepped, "output")
  prepped_data <- bind_cols(data.frame(classes) |> setNames(data |> select(all_of(class)) |> colnames()), 
                            map2(output, names(output), ~ .x |> set_names(paste(.y, names(.x), sep="."))))
  if(!is.null(residualised)) {
    prepped_data <- prepped_data |> bind_cols(data |> select(all_of(residualised)))
  }
  extra <- map(prepped, "extra")
  list(training = prepped_data,
       extra = extra)
}

# given variable information and a level map, do the mapping. Absent levels are mapped to 'new'
impute_score_ca0 <- function(var, extra) {
  var <- droplevels(var)
  var_levels <- pluck(extra, "var_levels")
  new.var_levels <- setdiff(levels(var), var_levels)
  X <- pluck(extra, "X")
  new_var_level_score <- X |> rownames_to_column("Var_Level") |> filter(Var_Level == "new") |> select(-Var_Level) |> as.list()
  test_score <- data.frame(Var_Level = var) |> left_join(X |> rownames_to_column("Var_Level"), by="Var_Level") |> replace_na(replace=new_var_level_score)
  list(test_score = test_score |> select(-Var_Level))
}

# iterate over the variable columns and use the extra info from before to remap the levels in the test data. 
prepare_test_ca0 <- function(data, extra, id, residualised=NULL) {
  # first remap the variable levels to the appropriate ordinal level
  id <- data |> select(all_of(id))
  var_cols <- data |> select(any_of(names(extra)))
  newdata_score <- map2(var_cols, extra, impute_score_ca0)
  output <- map(newdata_score,"test_score")
  newdata_pred <- map2_dfc(output, names(output), ~ .x |> set_names(paste(.y, names(.x), sep=".")))
  newdata_pred <- if(!is.null(residualised)) {
    bind_cols(id, data |> select(all_of(residualised)), newdata_pred)
  } else {
    bind_cols(id, newdata_pred)
  }
  newdata_pred
}
