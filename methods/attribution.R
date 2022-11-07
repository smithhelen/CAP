#### Source attribution of the SACNZ data (jejuni and coli) using the different methods ####

# Attribution 
attribution_fn <- function(Dat.train, Dat.test, d, method=ca0, ntrees=500, k=2, m=NULL, mp=100, axes=2){
  switch(method, 
         # prepare data using method of choice:
         ca0 = {
           train <- prepare_training_ca0(Dat.train, starts_with("CAMP"), "Source", axes=axes)
           test <- prepare_test_ca0(Dat.test, train$extra, "LabID")
         },
         pco = {
           train <- prepare_training_pco(Dat.train, starts_with("CAMP"), "Source", d, axes=axes)
           test <- prepare_test_pco(Dat.test, train$extra, "LabID")
         },
         cap = {  
           train <- prepare_training_cap(Dat.train, starts_with("CAMP"), "Source", d, k=k, m=m, mp=mp, axes=axes)
           test <- prepare_test_cap(Dat.test, train$extra, "LabID")
         }
  )       
  set.seed(3)
  # generate random forest models and make predictions for human data 
  rf_mod <- ranger(Source ~ ., data=train$training, oob.error = TRUE, num.trees=ntrees, respect.unordered.factors = TRUE)
  Prediction <- predict(rf_mod, data=test, predict.all = FALSE)$predictions
  out <- table(Prediction) %>% as.data.frame() # counts of predictions for each source
  out
}
