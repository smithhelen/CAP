#### functions to calculate misclassification rates of the individual trees according to the presence of absent levels 
#### as well as overall misclassification rates of the whole forests for the different methods ####

## NEED TO UPDATE FOR CAP

# functions
# (1) - calculate proportion correct predictions for trees with and without unique alleles
misclass_tree_fn <- function(dat){
  df <- dat |> 
    mutate(uses_unique = if_else(uses_unique == 0, "No", "Yes")) |>
    group_by(method,uses_unique,Source, prediction) |>
    summarise(n=n()) |>
    group_by(method,uses_unique,Source) |>
    mutate(N = sum(n)) |>
    mutate(prop = n/N)
  df
}

  
# (2) - calculate proportion correct classifications for forests (ie calculate misclassification rate)
MC_fn <- function(Dat.train, Dat.test, d=NULL, axes=2, mp=100, m=NULL, k=2, method=ca0, ntrees=500, residualised=NULL, id="LabID", class="Source"){
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
  uniques <- is_unique(Dat.test, train$extra) # can delete after testing
  set.seed(3)
  classes <- train$training |> pull({{class}})
  ranger_mod <- if(is.null({{residualised}})){
    ranger(classes ~ ., data=train$training |> select(-{{class}}), oob.error = TRUE, num.trees=ntrees, respect.unordered.factors = TRUE)
  } else {
    ranger(classes ~ ., data=train$training |> select(-{{class}}), oob.error = TRUE, num.trees=ntrees, respect.unordered.factors = "partition", 
           always.split.variables = {{residualised}})
  }
  preds_ranger <- predict(ranger_mod, data=test, predict.all = FALSE)$predictions
  df <- data.frame(id = test |> pull({{id}}), preds = preds_ranger) |> left_join(Dat.test |> rename(id = {{id}}) |> select(id, {{class}}))
  
  tree_preds_ranger <- predict(ranger_mod, data=test, predict.all = TRUE)$predictions |> as.data.frame() |> 
    mutate(id = test |> pull({{id}})) |> left_join(Dat.test |> rename(id = {{id}}) |> select(id, {{class}})) |> 
    pivot_longer(-c(id, {{class}}), names_to="tree", values_to="values") |> 
    left_join(data.frame(values = seq_along(ranger_mod$forest$levels), prediction = ranger_mod$forest$levels)) |>  # NOTE: ranger has a mod$forest$class.values that seems though it's not supposed to be used???
    mutate(tree = as.numeric(substring(tree, 2))) |> 
    select(id, tree, {{class}}, prediction)
  
  tree_preds_JM <- predict_by_tree(ranger_mod, test, uniques, id=id, residualised=residualised) |> 
    left_join(Dat.test |> rename(id = {{id}}) |> select(id, {{class}}))# can delete after testing
  
  out <- list(MC = df, tree_preds_ranger=tree_preds_ranger, tree_preds_JM = tree_preds_JM)# replace with 'df' after testing
  out
}

calc_misclassification <- function(df) {
  # weights for each fold
  w <- df |> group_by(Fold) |> tally() |> pull(n)
  
  # misclassification rates of the folds
  lm <- df |> 
    mutate(wrong = truths != preds) |>
    group_by(Fold) |>
    summarise(miss = sum(wrong)/n()) |> 
    lm(miss ~ 1, weights = w, data=_) |> summary()
  
  # weighted mean and standard error of the misclassification rates
  av <- coef(lm)[,"Estimate"]
  se <- coef(lm)[,"Std. Error"]
  
  # confusion matrix
  conf <- df |> 
    group_by(Fold, truths, preds) |> 
    summarise(n=n()) |> 
    mutate(N=sum(n), p=n/sum(n)) |> 
    ungroup() |>
    mutate(w = as.numeric(paste(factor(Fold, labels=w))), t_p = paste(truths, preds, sep="_")) |> 
    split(~t_p)

  # weighted mean and standard error of the confusion matrices
  conf.lm <- map(conf, function(x) {lm(p~1, weights=w, data=x) |> summary()})
  conf.av <- map(conf.lm, function(x) coef(x)[,"Estimate"])
  conf.se <- map(conf.lm, function(x) coef(x)[,"Std. Error"])

  out <- list(av=av, se=se, conf.av=conf.av, conf.se=conf.se)
  out
}

misclass_fn <- function(Dat.train, Dat.test, method=ca0, class="Source", id="LabID", d=NULL, k=2, m=NULL, mp=100, axes=2, ntrees=500, residualised=NULL){
  DF <- map2_dfr(Dat.train, Dat.test, ~MC_fn(.x,.y,method={{method}},d=d, k=k, m=m, mp=mp, axes=axes, ntrees=ntrees, residualised=residualised), .id="Fold")
  DF# can delete after testing
  #answer <- calc_misclassification(DF) #change back after testing
  #answer
}
