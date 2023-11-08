# Pull out individual tree predictions

is_unique_level <- function(var, extra) {
  if(!is.numeric(var))  {var <- droplevels(var)}
  var_levels <- pluck(extra, "var_levels")
  is_unique <- !(var %in% var_levels)
  dim <- pluck(extra, "dim")
  suffix <- pluck(extra, "suffix")
  out <- data.frame(matrix(is_unique, ncol=dim, nrow=length(is_unique))) %>% set_names(suffix)
  out
}

is_unique <- function(data, list_of_extras) {
  var_cols <- data %>% select(any_of(names(list_of_extras)))
  output <- map2(var_cols, list_of_extras, is_unique_level)
  uniques <- map2_dfc(output, names(output), 
                      ~ if(!is.null(names(.x))) {.x %>% set_names(paste(.y, names(.x), sep="."))} 
                      else(.x %>% set_names(paste(.y))))
  uniques
}

# do the tree prediction
predict_row <- function(tree, data_row, uniques_row) {
  # pass the data down the tree row by row
  names(uniques_row) <- names(uniques_row) %>% sub("\\..*","",.)
  uses_unique <- 0
  prediction <- NULL
  vars_used_in_tree <- NULL
  unique_vars_used_in_tree <- NULL
  row <- 1 # This code assumes that nodeID = row-1
  while (TRUE) {
    # split data to the left and right nodes
    if (tree$terminal[row]) {
      # terminal node - we have our answer
      prediction = tree$prediction[row]
      break;
    }
    # is our level unique?
    split = tree$splitvarName[row]  #name of var used in tree
    # check if the named split variable is in the uniques row
    name_of_split_variable <- split %>% sub("\\..*","",.)
    if (name_of_split_variable %in% names(uniques_row) && #' This variable is in the uniques row
        uniques_row[[name_of_split_variable]]) {
      uses_unique = uses_unique + 1
      unique_vars_used_in_tree <- c(unique_vars_used_in_tree, name_of_split_variable)
    }
    # go down the tree
    if (data_row[[split]] <= tree$splitval[row]) { # Ranger uses <= here, and this gives same result
      # left tree
      row <- tree$leftChild[row] + 1
    } else {
      # right tree
      row <- tree$rightChild[row] + 1
    }
    vars_used_in_tree <- c(vars_used_in_tree, name_of_split_variable)
  }
  tibble(prediction=prediction, uses_unique=uses_unique, splitting_vars=list(vars_used_in_tree), 
         unique_splitting_vars=list(unique_vars_used_in_tree))
}

my_treeInfo <- function(mod, tree_number) {
  tree <- treeInfo(mod, tree_number)
  # treeInfo produces a splits column with both numeric and categorical splits, where
  # the categorical splits are decoded from their bitpattern encoding for readability.
  # We want the bitpattern encoding as that is the most efficient way of matching against
  # splits, so we'll override it here:
  tree$splitval = mod$forest$split.values[[tree_number]]
  tree
}

predict_tree <- function(mod, tree_number, nd, nu, id) {
  tree <- my_treeInfo(mod, tree_number)
  out_dfr <- map2_dfr(nd, nu, ~predict_row(tree, .x, .y))
  out_dfr |>
    mutate(tree = tree_number,
           id = id)
}

# do the predictions
predict_by_tree <- function(mod, new_data, new_unique, id) {
  nd <- split(new_data, 1:nrow(new_data))  # list, each entry is a row of test data
  nu <- split(new_unique, 1:nrow(new_unique))  # list, each entry is a row of uniques (ie TRUE or FALSE for each var)
  id = new_data %>% pull(id)
  predictions <- map_dfr(seq_len(mod$num.trees), ~predict_tree(mod=mod, tree_number=., nd=nd, nu=nu, id=id))
  predictions
}
