#### Prepare test data and training data for ca method ####

# Function to map from variable levels to scores
factor_to_ca_score <- function(var, class, k) {
  var_levels <- droplevels(var)
  if (nlevels(var_levels) < 2) {
    # if we only have one var_level we can't do anything
    return(list(output = rep(1, times=length(var)), 
                extra = list(var_levels=levels(var_levels), dim = 1, suffix=NULL, 
                             score = data.frame(Var_Level = levels(var_levels), Rank = 1))))
  }
  ct <- table(Var_Level=var_levels, class=class)
  if(is.null(k)){k <- ncol(ct)-1}
  # add zero row - can't add zeros to ct as P will complain. So add 1/num.classes to P (equal probabilities across classes)
  new <- matrix(rep((1/nlevels(class)),times=nlevels(class)),nrow=1, 
                ncol=nlevels(class), dimnames = list(Var_Level="new",Class=colnames(ct)))
  P <- rbind(ct/rowSums(ct),new)
  S <- cov.wt(P, wt = c(rowSums(ct),0))$cov
  
  eigen_S <- eigen_decomp(S, symmetric=TRUE) ## PCA of weighted covariance matrix of class probabilites
  # Restrict to a maximum of eigenvectors set by "k" (the number of axes) (default is NULL = ncol(ct)-1)
  nlambdas <- min(sum(eigen_S$values > epsilon), k)
  pc <- eigen_S$vectors
  X <- P %*% pc[, seq_len(nlambdas), drop=FALSE] |> as.data.frame() 
  score <- left_join(data.frame(Var_Level = var_levels), X |> rownames_to_column("Var_Level"), by = "Var_Level")
  output <- score |> select(-Var_level)
  list(output = output,
       extra = list(X=X, k=k, var_levels=levels(var_levels), dim = ncol(X), suffix = colnames(X)))
}

# iterate over the variable columns and grab the output as a new data.frame for ca, and store the absent level stuff for later 
prepare_training_ca <- function(data, var_cols, class, k=NULL) {
  # pull out our variable columns and class
  var_cols <- dplyr::select(data, all_of(var_cols))
  classes   <- data |> pull(all_of(class))
  # iterate over the variable columns, and convert
  prepped <- map(var_cols, factor_to_ca_score, class = classes, k=k)
  output <- map(prepped,"output")

  prepped_data <- bind_cols(data.frame(classes) |> setNames(data |> select(all_of(class)) |> colnames()), 
                            map2(output, names(output), ~ .x |> set_names(paste(.y, names(.x), sep="."))))
  extra <- map(prepped, "extra")
  list(training = prepped_data,
       extra=extra)
}

# score absent levels as infinite
impute_ordinal_ca <- function(var, extra) {
  var <- droplevels(var)
  var_levels <- pluck(extra, "var_levels")
  new.var_levels <- setdiff(levels(var), var_levels)
  X <- pluck(extra, "X")

  
}

# given var information and a level map, do the mapping. New var_levels are mapped to 'new'
impute_score_ca <- function(var, extra) {
  X <- pluck(extra, "X")
  new_var_level_score <- X %>% filter(Var_level == "new") %>% pull(PC1)
  test_score <- data.frame(Var_level = var) %>%
    left_join(X, by="Var_level") %>%
    replace_na(list(PC1 = new_var_level_score))
  test_score %>% pull(PC1)
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
prepare_test_ca <- function(data, list_of_extras, id) {
  # first remap the variable levels to the appropriate ordinal level
  id <- data |> pull({{id}})
  test_data <- data |> select(any_of(names(list_of_extras)))
  newdata_pred <- map2_dfc(test_data, list_of_extras, impute_score_ca)
  newdata_pred <- bind_cols(id=id, newdata_pred)
  newdata_pred
}
