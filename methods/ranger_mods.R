# function to prepare the training and testing data
prep_data <- function(Dat.train, Dat.test, method="ca0", residualised=NULL, d=NULL, id, class, var_id, k=NULL, m=NULL, mp=99, c=NULL){
  switch(method,
         ca = {
           train <- prepare_training_ca(Dat.train, starts_with(var_id), class=class)
           test <- prepare_test_ca(Dat.test, train$extra, id=id)
         },
         binary = {
           train <- prepare_training_binary(Dat.train, starts_with(var_id), class=class)
           test <- prepare_test_binary(Dat.test, train$extra, id=id)
         },
         ca0 = {
           train <- prepare_training_ca0(Dat.train, starts_with(var_id), class=class, k=k, residualised=NULL)
           test <- prepare_test_ca0(Dat.test, train$extra, id=id)
         },
         pco = {
           train <- prepare_training_pco(Dat.train, starts_with(var_id), class=class, d=d, m=m, mp=mp, residualised=NULL)
           test <- prepare_test_pco(Dat.test, train$extra, id=id)
         },
         cap = {
           train <- prepare_training_cap(Dat.train, starts_with(var_id), class=class, d=d, k=k, m=m, mp=mp, c=c, residualised=NULL)
           test <- prepare_test_cap(Dat.test, train$extra, id=id)
         },
         similarity = {
           train <- prepare_training_similarity(Dat.train, starts_with(var_id), class=class, d=d)
           test <- prepare_test_similarity(Dat.test, train$extra, id=id)
         }
  )
  prepped <- list(train=train, test=test)
}

# function to prepare the training and testing data and train a random forest (ranger) model
prep_and_train <- function(Dat.train, Dat.test, method="ca0", residualised=NULL, d=NULL, id, class, var_id, k=NULL, m=NULL, mp=99, c=NULL, ntrees=500){
  # prepare data
  prepped <- prep_data(Dat.train, Dat.test, method=method, residualised=residualised, d=d, id=id, class=class, var_id=var_id, k=k, m=m, mp=mp, c=c)
  # train model
  classes <- prepped$train$training |> pull(all_of({{class}}))
  ranger_mod <- if(method=="binary") {
    ranger(dependent.variable.name = class, data = prepped$train$training, oob.error = TRUE, num.trees=ntrees, respect.unordered.factors = TRUE)
  } else {
    ranger(classes ~ ., data=prepped$train$training |> select(-any_of({{class}})), oob.error = TRUE, num.trees=ntrees, respect.unordered.factors = TRUE)
  }
  model.info <- list(train=prepped$train, test=prepped$test, class=class, ranger_mod=ranger_mod, id={{id}})
  model.info
}
