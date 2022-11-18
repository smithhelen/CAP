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
predict_row <- function(tree, data_row, uniques_row, residualised) {
  # pass the data down the tree row by row
  uses_unique <- 0
  prediction <- NULL
  vars_used_in_tree <- NULL
  unique_vars_used_in_tree <- NULL
  uses_res_var <- 0
  row <- 1 # This code assumes that nodeID = row-1
  while (TRUE) {
    #    cat("row ", row, "\n")
    # split data to the left and right nodes
    if (tree$terminal[row]) {
      # terminal node - we have our answer
      prediction = tree$prediction[row]
      break;
    }
    # is our level unique?
    split = tree$splitvarName[row]  #name of var used in tree
    # cat(split)   # for checking errors
    if(!is.null(residualised) && split == residualised){
      # NOTE: This code assumes residualised is a factor. It doesn't have to be, but
      #       if it's not then this will break. We could probably fix this by appending another
      #       column to tree with whether the split is an unordered variable or not. If it's
      #       ordered it's a factor with gold-standard (bitwise matching) so should use
      #       this branch. Else, we can use the other one.
      # OK, let's figure out if we hit it. We hit it if 2^(split-1) is in our splitval
      uses_res_var = uses_res_var + 1
      # go down the tree. we use bitwise And here as the right branch is encoded as a bit pattern
      if (bitwAnd(data_row[[split]], tree$splitval[row])) {
        # right tree
        row <- tree$rightChild[row] + 1
      } else {
        # left tree
        row <- tree$leftChild[row] + 1
      }
    } else {
      if (uniques_row[[split]]) {
        uses_unique = uses_unique + 1
        unique_vars_used_in_tree <- c(unique_vars_used_in_tree, split |> sub(pattern = "\\..*", replacement = ""))
      }
      # go down the tree
      if (data_row[[split]] <= tree$splitval[row]) { # Ranger uses <= here, and this gives same result
        # left tree
        row <- tree$leftChild[row] + 1
      } else {
        # right tree
        row <- tree$rightChild[row] + 1
      }
    }
    vars_used_in_tree <- c(vars_used_in_tree, split |> sub(pattern = "\\..*", replacement = ""))
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

predict_tree <- function(mod, tree_number, nd, nu, id, residualised) {
  #cat("working on tree", tree_number, "\n")
  tree <- my_treeInfo(mod, tree_number)
  out_dfr <- map2_dfr(nd, nu, ~predict_row(tree, .x, .y, residualised=residualised))
  out_dfr |>
    mutate(tree = tree_number,
           id = id)
}

# do the predictions
predict_by_tree <- function(mod, new_data, new_unique, id, residualised) {
  x <- new_data
  if (!is.null(residualised)) {
    # need to replace our (factor) variable residualised with the corresponding
    # integers as per ranger. Ranger uses data.matrix for this (see ranger:::predict.ranger.forest)
    # and then encodes the integers as 2^(number-1) so it can bit-wise match
    # the bitpattern for 'going right' down the tree. We reproduce that here.
    
    # Technically this doesn't need to happen in the case of residualised not being a factor though.
    # We don't cover that case here yet, but we probably could by using mod$is.ordered: if that is false
    # for the 'residualised' entry then we need to do this magic. If it is true, we don't, and can
    # tree CC the same as always.
    x[[residualised]] <- as.integer(x[[residualised]])
    x[[residualised]] <- 2^(x[[residualised]]-1)
  }
  
  nd <- split(x, 1:nrow(x))  # list, each entry is a row of test data
  nu <- split(new_unique, 1:nrow(new_unique))  # list, each entry is a row of uniques (ie TRUE or FALSE for each var)
  id = new_data |> pull({{id}})
  predictions <- map_dfr(seq_len(mod$num.trees), ~predict_tree(mod=mod, tree_number=., nd=nd, nu=nu, id=id, residualised=residualised))
  
  # BUG IN RANGER. treeInfo() produces incorrect forest levels 
  if (!is.null(mod$forest$levels)) {
    # fix up the levels. We have to undo the command:
    # factor(result$prediction, levels = forest$class.values, labels = forest$levels)
    # there is the following new code in ranger GitHub but the package is not yet updated (10/11/2022)
    #if (!is.null(forest$levels)) { result$prediction <- integer.to.factor(result$prediction, labels = forest$levels)  }
    fixup <- data.frame(predict = mod$forest$levels[mod$forest$class.values], treeinfo = mod$forest$levels)
    predictions <- predictions |>
      left_join(fixup, by=c("prediction" = "treeinfo")) |>
      select(id, tree, prediction = predict, uses_unique, splitting_vars, unique_splitting_vars)
  }
  predictions
}

# pull out individual tree decisions for each observation
tree_fn <- function(Dat.train, Dat.test, d=NULL, axes=2, mp=100, m=NULL, k=2, method=ca0, ntrees=500, residualised=NULL, id="LabID", class="Source"){
  switch(method, 
         ca0 = {
           train <- prepare_training_ca0(Dat.train, starts_with("CAMP"), class="Source", axes=axes, residualised=residualised)
           test <- prepare_test_ca0(Dat.test, train$extra, id={{id}}, residualised=residualised)
         },
         pco = {
           train <- prepare_training_pco(Dat.train, starts_with("CAMP"), class="Source", d, axes=axes, residualised=residualised)
           test <- prepare_test_pco(Dat.test, train$extra, id={{id}}, residualised=residualised)
         },
         cap = {  
           train <- prepare_training_cap(Dat.train, starts_with("CAMP"), class="Source", d, axes=axes, k=k, m=m, mp=mp, residualised=residualised)
           test <- prepare_test_cap(Dat.test, train$extra, id={{id}}, residualised=residualised)
         }
  )       
  uniques <- is_unique(Dat.test, train$extra)
  set.seed(3)
  classes <- train$training |> pull({{class}})
  rf_mod <- if(is.null(residualised)){
    ranger(classes ~ ., data=train$training |> select(-{{class}}), oob.error = TRUE, num.trees=ntrees, respect.unordered.factors = TRUE)
  } else {
    ranger(classes ~ ., data=train$training |> select(-{{class}}), oob.error = TRUE, num.trees=ntrees, respect.unordered.factors = "partition", 
           always.split.variables = residualised)
  }
  #list(rf_mod=rf_mod, test=test, uniques=uniques, id=id) # for troubleshooting
  tree_preds <- predict_by_tree(rf_mod, test, uniques, id=id, residualised=residualised)
  answer <- tree_preds  |> left_join(Dat.test |> rename(id = {{id}}) |> select(id, {{class}}))
  answer
}



## make predictions for Human data
#predictions <- predict(rf_mod, data=test, predict.all = FALSE)$predictions
#table(predictions) %>% as.data.frame() # counts of predictions for each source



