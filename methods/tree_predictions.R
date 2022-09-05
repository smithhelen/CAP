# Pull out individual tree predictions

is_unique_level <- function(var, extra) {
  var_levels <- pluck(extra, "var_levels")
  var <- droplevels(var)
  is_unique <- !(var %in% var_levels)
  num_vars <- pluck(extra, "num_vars")
  var_names <- pluck(extra, "var_names")
  out <- data.frame(matrix(is_unique, ncol=num_vars, nrow=length(is_unique))) |> set_names(var_names)
  out
}

is_unique <- function(data, extra) {
  var_cols <- data |> select(any_of(names(extra)))
  output <- map2(var_cols, extra, is_unique_level)
  uniques <- map2_dfc(output, names(output), 
                      ~ if(!is.null(names(.x))) {.x |> set_names(paste(.y, names(.x), sep="."))} 
                      else(.x |> set_names(paste(.y))))
  uniques
}

# do the tree prediction
predict_row <- function(tree, data_row, uniques_row) {
  # pass the data down the tree row by row
  names(uniques_row) <- names(uniques_row) |> sub("\\..*","",.)
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
    #cat(split)   # for checking errors
    if (uniques_row[[split |> sub("\\..*","",.)]]) {
      uses_unique = uses_unique + 1
      unique_vars_used_in_tree <- c(unique_vars_used_in_tree, split |> sub("\\..*","",.))
    }
    # go down the tree
    if (data_row[[split]] <= tree$splitval[row]) { # Ranger uses <= here, and this gives same result
      # left tree
      row <- tree$leftChild[row] + 1
    } else {
      # right tree
      row <- tree$rightChild[row] + 1
    }
    vars_used_in_tree <- c(vars_used_in_tree, split |> sub("\\..*","",.))
  }
  tibble(prediction=prediction, uses_unique=uses_unique, splitting_vars=list(vars_used_in_tree), unique_splitting_vars=list(unique_vars_used_in_tree))
}

predict_tree <- function(mod, tree_number, nd, nu, id) {
  #cat("working on tree", tree_number, "\n")
  tree <- treeInfo(mod, tree_number)
  out_dfr <- map2_dfr(nd, nu, ~predict_row(tree, .x, .y))
  out_dfr |>
    mutate(tree = tree_number,
           row = id)
}

# do the predictions
predict_by_tree <- function(mod, new_data, new_unique) {
  nd <- split(new_data, 1:nrow(new_data))  # list, each entry is a row of test data
  nu <- split(new_unique, 1:nrow(new_unique))  # list, each entry is a row of uniques (ie TRUE or FALSE for each var)
  id = new_data |> pull(id)
  predictions <- map_dfr(seq_len(mod$num.trees), ~predict_tree(mod=mod, tree_number=., nd=nd, nu=nu, id=id))
  
  # BUG IN RANGER. treeInfo() produces incorrect forest levels 
  if (!is.null(mod$forest$levels)) {
    # fix up the levels. We have to undo the command:
    # factor(result$prediction, levels = forest$class.values, labels = forest$levels)
    fixup <- data.frame(predict = mod$forest$levels[mod$forest$class.values], treeinfo = mod$forest$levels)
    predictions <- predictions |>
      left_join(fixup, by=c("prediction" = "treeinfo")) |>
      select(row, tree, prediction = predict, uses_unique, splitting_vars, unique_splitting_vars)
  }
  predictions
}

