#### Prepare test data and training data for ca method ####

# Function to map from variable levels to scores
factor_to_ordinal_ca <- function(var, class) {
  var_levels <- droplevels(var)
  if (nlevels(var_levels) < 2) {
    # if we only have one var_level we can't do anything
    return(list(output = rep(1, times=length(var)), 
                extra = list(var_levels=levels(var_levels), dim = 1, suffix=NULL, 
                             score = data.frame(Var_Level = levels(var_levels), Rank = 1))))
  }
  ct <- table(Var_Level=var_levels, class=class)
  P <- ct/rowSums(ct)
  S <- cov.wt(P, wt = rowSums(ct))$cov
  pc1 <- eigen(S)$vectors[,1]   ## PCA of weighted covariance matrix of class probabilites
  score <- P %*% pc1 %>% as.data.frame() %>% rownames_to_column("Var_Level") %>% setNames(., c("Var_Level","PC1")) %>% mutate(Rank = rank(PC1, ties.method = "first"))
  output_factor <- data.frame(Var_Level = var_levels) %>% left_join(score, by="Var_Level") %>% pull(Rank)
  list(output = output_factor, 
       extra = list(score = score, var_levels = levels(var_levels), dim = 1, suffix=NULL))
}

# iterate over the variable columns and grab the output as a new data.frame for ca, and store the absent level stuff for later 
prepare_training_ca <- function(data, var_cols, class) {
  # pull out our variable columns and class
  var_cols <- dplyr::select(data, all_of(var_cols))
  classes   <- data %>% pull(class)
  # iterate over the variable columns, and convert
  prepped <- map(var_cols, factor_to_ordinal_ca, class = classes)
  prepped_data <- bind_cols(data.frame(classes) %>% setNames(data %>% select(all_of(class)) %>% colnames), 
                            map_dfc(prepped, "output"))
  extra <- map(prepped, "extra")
  list(training = prepped_data,
       extra=extra)
}

# score absent levels as infinite
impute_ordinal_ca <- function(var, extra) {
  score <- pluck(extra, "score")
  new_var_level_rank <- score %>% pull(Rank) %>% max() + 1
  ranked_var <- data.frame(Var_Level = var) %>%
    left_join(score, by="Var_Level") %>%
    replace_na(list(Rank = new_var_level_rank))
  ranked_var %>% pull(Rank)
}

# iterate over the variable columns and use the extra info from before to remap the levels in the test data. 
prepare_test_ca <- function(data, list_of_extras, id) {
  # first remap the variable levels to the appropriate ordinal level
  id <- data %>% pull({{id}})
  test_data <- data %>% select(any_of(names(list_of_extras)))
  newdata_pred <- map2_dfc(test_data, list_of_extras, impute_ordinal_ca)
  newdata_pred <- bind_cols(id=id, newdata_pred)
  newdata_pred
}
