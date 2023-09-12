#### Prepare test data and training data for binary method ####

# load libraries
library(varhandle) # for to.dummy()

# iterate over the variable columns and grab the output as a new data.frame, and store the absent level stuff for later 
binary_map <- function(var, var_name) {
  var <- droplevels(var)
  var_binary <- data.frame(to.dummy(var, var_name))
  list(output = var_binary, 
       extra = list(var_levels=levels(var), dim=1, var_name=var_name, suffix=NULL))
}

prepare_training_binary <- function(data, var_cols, class) {
  # pull out the variable columns and class
  var_cols <- dplyr::select(data, all_of(var_cols))
  classes   <- data %>% pull(class)
  # iterate over the variable columns, and convert
  prepped <- map2(var_cols, names(var_cols), binary_map)
  prepped_data <- bind_cols(data.frame(classes) %>% setNames(data %>% select(all_of(class)) %>% colnames), 
                            map_dfc(prepped, "output"))
  extra <- map(prepped, "extra")
  list(training = prepped_data,
       extra=extra)
}

# score absent levels as zero
impute_binary <- function(var, extra) {
  var_levels <- pluck(extra, "var_levels")
  var_name <- pluck(extra, "var_name")
  vars <- paste(var_name, var_levels, sep=".")
  var_binary <- data.frame(to.dummy(var, var_name)) %>% select(any_of(vars))
  missing_vars <- setdiff(vars, colnames(var_binary))
  empty_cols <- data.frame(matrix(0, nrow=nrow(var_binary),ncol=length(missing_vars)))
  colnames(empty_cols) <- missing_vars
  test_binary <- var_binary %>% bind_cols(empty_cols) %>% select(any_of(vars))
  # sanity check that the columns we have are only those that appear in the training data (and that
  # we have all of them in the right order)
  stopifnot(names(test_binary) == vars)
  return(test_binary)
}

# iterate over the variable columns and use the extra info from before to remap the levels in the test data
prepare_test_binary <- function(data, list_of_extras, id) {
  test_data <- data %>% select(any_of(names(list_of_extras)))
  newdata_pred <- map2_dfc(test_data, list_of_extras, impute_binary)
  newdata_pred <- bind_cols(data |> select(all_of(id)), newdata_pred)
  newdata_pred
}






















