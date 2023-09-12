#' functions to calculate misclassification rates of the individual trees according to the presence of absent levels 
#' as well as overall misclassification rates of the whole forests for the different methods 

source("methods/tree_predictions.R")

# (1) - calculate proportion correct predictions for trees with and without unique alleles
tree_predictions <- function(data, model){
  prep.train <- pluck(model, "train")
  prep.test <- pluck(model, "test")
  class <- pluck(model, "class")
  id <- pluck(model, "id")
  ranger_mod <- pluck(model, "ranger_mod")
  uniques <- is_unique(data, prep.train$extra)
  predictions <- predict_by_tree(ranger_mod, new_data=prep.test, new_unique=uniques, id={{id}})
  tree_preds <- predictions |> left_join(data |> rename(id = {{id}}) |> select(id, any_of({{class}})), by="id")
  tree_preds
}

tree_misclassifications <- function(testdata, model, class){
  tree_preds <- map2_df(testdata, model, tree_predictions, .id = "fold")
 
  # misclassification rates of the folds
  tree.res <- tree_preds |> 
    mutate(uses_unique = if_else(uses_unique == 0, "No", "Yes")) |>
    group_by(fold, uses_unique, across(class), prediction) |>
    summarise(n=n()) |>
    group_by(fold, uses_unique, across(class)) |>
    mutate(N = sum(n)) |>
    mutate(prop = n/N)
  
  # overall misclassification rate
  tree.mc <- tree_preds |> mutate(uses_unique = if_else(uses_unique == 0, "No", "Yes")) |>
    group_by(uses_unique, across(class), prediction) |>
    summarise(n=n()) |>
    group_by(uses_unique, across(class)) |>
    mutate(N = sum(n)) |>
    mutate(prop = n/N)
  tree.yes <- tree.mc |> ungroup() |> filter(uses_unique=="Yes") |> arrange(desc(prop)) |> 
    select(all_of(class), prediction, prop) 
  tree.no <- tree.mc |> ungroup() |> filter(uses_unique=="No") |> arrange(desc(prop)) |> 
    select(all_of(class), prediction, prop)
  tree.mc <- right_join(tree.yes, tree.no, by=c({{class}}, "prediction"), suffix=c(".yes",".no"))
  
  out=list(tree.res=tree.res, tree.mc=tree.mc)
  out
}


# (2) - calculate proportion correct classifications for forests (ie calculate misclassification rate)
forest_predictions <- function(data, model){
  prep.test <- pluck(model, "test")
  class <- pluck(model, "class")
  id <- pluck(model, "id")
  ranger_mod <- pluck(model, "ranger_mod")
  predictions <- predict(ranger_mod, data=prep.test, predict.all = FALSE)$predictions
  forest_preds <- data.frame(id = prep.test |> pull(id), prediction = predictions) |> 
    left_join(data |> rename(id = {{id}}) |> select(id, any_of({{class}})), by="id")
  forest_preds
}


forest_misclassifications <- function(testdata, model, class) {
  # forest predictions
  forest_preds <- map2_df(testdata, model, forest_predictions, .id = "fold")
  
  # weights for each fold
  w <- forest_preds |> group_by(fold) |> tally() |> pull(n)
  
  # misclassification rates of the folds
  lm <- forest_preds |> 
    rename(truths = {{class}}) |> 
    mutate(wrong = truths != prediction) |>
    group_by(fold) |>
    summarise(miss = sum(wrong)/n()) |> 
    lm(miss ~ 1, weights = w, data=_) |> summary()
  
  # weighted mean and standard error of the misclassification rates
  av <- coef(lm)[,"Estimate"]
  se <- coef(lm)[,"Std. Error"]
  
  # confusion matrix
  conf <- forest_preds |> 
    rename(truths = {{class}}) |> 
    group_by(fold, truths, prediction) |> 
    summarise(n=n()) |> 
    mutate(N=sum(n), p=n/sum(n)) |> 
    ungroup() |>
    mutate(w = as.numeric(paste(factor(fold, labels=w))), t_p = paste(truths, prediction, sep="_")) |> 
    split(~t_p)
  
  # weighted mean and standard error of the confusion matrices
  conf.lm <- map(conf, function(x) {lm(p~1, weights=w, data=x) |> summary()})
  conf.av <- map(conf.lm, function(x) coef(x)[,"Estimate"])
  conf.se <- map(conf.lm, function(x) coef(x)[,"Std. Error"])
  
  # overall misclassification rate
  forest.mc <- forest_preds |> group_by(across(class), prediction) |> 
    summarise(n=n()) |> mutate(N = sum(n)) |> mutate(prop = n/N) |> select(-n,-N) 
  
  out <- list(av=av, se=se, conf.av=conf.av, conf.se=conf.se, forest.mc=forest.mc)
  out
}

mc_function <- function(Dat.test.list, model.info.list, class, tree=TRUE, forest=TRUE){
  if(forest){
    # forest predictions
    forest.mc <- forest_misclassifications(Dat.test.list, model.info.list, class={{class}})$forest.mc
    results <- forest.mc
  }
  
  if(tree){
    # tree predictions
    tree.mc <- tree_misclassifications(Dat.test.list, model.info.list, class={{class}})$tree.mc
    results <- tree.mc
  }
  
  if(forest & tree){
    # merge forest and tree results
    results <- left_join(tree.mc, forest.mc, by=c(class, "prediction")) |> 
      rename(tree.absent=prop.yes,tree.no.absent=prop.no,forest=prop)
  }
  
  # output
  results
}
