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
MC_fn <- function(Dat.train, Dat.test, d, axes=1, mp=99, method=ca0, ntrees=500){
  switch(method, 
         ca0 = {
           train <- prepare_training_ca0(Dat.train, starts_with("CAMP"), "Source", axes=axes)
           test <- prepare_test_ca0(Dat.test, train$extra, "LabID")
         },
         pco = {
           train <- prepare_training_pco(Dat.train, starts_with("CAMP"), "Source", d, axes=axes)
           test <- prepare_test_pco(Dat.test, train$extra, "LabID")
         },
         cap = {  
           train <- prepare_training_cap(Dat.train, starts_with("CAMP"), "Source", d, axes=NULL, mp=mp)
           test <- prepare_test_cap(Dat.test, train$extra, "LabID")
         }
  )       
  set.seed(3)
  ranger_mod <- ranger(Source ~ ., data=train$training, oob.error = TRUE, num.trees=ntrees, respect.unordered.factors = TRUE)
  pred <- predict(ranger_mod, data=test, predict.all = FALSE)$predictions
  df <- data.frame(preds = pred, truths = Dat.test |> pull(Source))
  df
}

calc_misclassification <- function(df) {
  # weights for each fold
  w <- df |> group_by(Fold) |> tally() |> pull(n)
  
  # misclassification rates of the folds
  lm <- df |> 
    mutate(wrong = truths != preds) |>
    group_by(Fold) |>
    summarise(miss = sum(wrong)/n()) |> 
    lm(miss ~ 1, weights = w, data=.) |> summary()
  
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

misclass_fn <- function(Dat.train, Dat.test, method=ca0, d=NULL, axes=2){
  switch(method,
         cap = {
           cat("Not ready for CAP yet")
           return()
         },
         {
           DF <- map2_dfr(Dat.train, Dat.test, ~MC_fn(.x,.y,method={{method}},d=d, axes=axes), .id="Fold")
           answer <- calc_misclassification(DF)
           answer
         }
  )
}
