#### Prepare test data and training data for Cerda's similarity encoding method ####

# Function to map from variable levels to scores.
# Output is both the score information for the variable (i.e. a vector of the same length as the input vector) and
# the level mapping data.frame (mapping from var_level to score)
factor_to_similarity_score <- function(var, dist) {
  var_levels <- droplevels(var)
  d.train <- dist[levels(var_levels),levels(var_levels),drop=FALSE]
  score <- left_join(data.frame(Var_Level = var_levels), data.frame(d.train) |> rownames_to_column("Var_Level"), by = "Var_Level")
  output <- score |> select(-Var_Level)
  list(output = output,
       extra = list(d=dist, var_levels=levels(var_levels), suffix = NULL, dim=1))
}

prepare_training_similarity <- function(data, var_cols, class, d) {
  # pull out our var_cols and class
  var_cols <- dplyr::select(data, all_of(var_cols))
  classes   <- data %>% pull(class)
  # iterate over the var columns and distance matrices, and convert
  prepped <- map2(var_cols, d, factor_to_similarity_score) %>% compact() # removes empties
  output <- map(prepped,"output")
  prepped_data <- bind_cols(data.frame(classes) %>% setNames(data %>% select(all_of(class)) %>% colnames), 
                            map2(output, names(output), ~ .x %>% set_names(paste(.y, names(.x), sep="."))))
  extra <- map(prepped, "extra")
  list(training = prepped_data,
       extra = extra)
}


# given var information and extra information, do the mapping. 
impute_score_similarity <- function(var, extra) {
  var <- droplevels(var)
  var_levels <- pluck(extra, "var_levels")
  new.var_levels <- setdiff(levels(var), var_levels)
  d <- pluck(extra, "d")
  new_scores <- d[c(var_levels, new.var_levels), var_levels]
  var_level_score <- new_scores %>% as.data.frame() %>% rownames_to_column("Var_Level")
  output <- data.frame(Var_Level = var) %>% left_join(var_level_score, by = "Var_Level") %>% select(-Var_Level)
  list(output = output)
}

prepare_test_similarity <- function(data, list_of_extras, id) {
  test_data <- data %>% select(any_of(names(list_of_extras)))
  prepped <- map2(test_data, list_of_extras, impute_score_similarity)
  output <- map(prepped, "output")
  newdata_pred <- map2_dfc(output, names(output), ~ .x %>% set_names(paste(.y, names(.x), sep=".")))
  newdata_pred <- bind_cols(data |> select(all_of(id)), newdata_pred)
  newdata_pred
}

