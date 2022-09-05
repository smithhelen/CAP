### pull out variable importance values ###
# note this is not the best way to determine the importance of variables 
# as it uses the OOB data but is a simple way to check if the same data is being used for the different methods (not CAP yet)
# use on whole training at once (as per attribution function)

var_impt_fn <- function(Dat.train, d=NULL, axes=2, mp=99, method=ca0, ntrees=500){
  switch(method, 
         ca0 = {
           train <- prepare_training_ca0(Dat.train, starts_with("CAMP"), "Source", axes=axes)
         },
         pco = {
           train <- prepare_training_pco(Dat.train, starts_with("CAMP"), "Source", d, axes=axes)
         },
         cap = {
           train <- prepare_training_cap(Dat.train, starts_with("CAMP"), "Source", d, mp=mp, axes=axes)
         }
  )       
  set.seed(3)
  ranger_mod <- ranger(Source ~ ., data=train$training, oob.error = TRUE, num.trees=ntrees, respect.unordered.factors = TRUE, importance="permutation")
  var_impt <- ranger_mod$variable.importance
  top10 <- var_impt |> as.data.frame() |> rownames_to_column("gene") |> rename(var_impt = ".") |> arrange(desc(var_impt))
  top10
}

