# Load libraries and functions
source("methods/libs_fns.R")
source("methods/tree_predictions.R")


# Load sample data 
# data is genetic data with alleles (levels) from 10 genes (variables) across 3 sources (classes), n=10 for each source
load("../CAP_data/data/sample_run_data.RData") # sample data
load("../CAP_data/data/sample_run_SeqDat.RData") # sequence information for sample data



# Prepare data
# calculate hamming distances between levels of predictor variable
source("methods/hamming.R")
genes <- names(SeqDat |> select(starts_with("Var")))
list_of_distance_matrices <- map(genes, ~ dfun(gene=., dat=LostInTheForest, seqdat=SeqDat))
names(list_of_distance_matrices) <- genes
list_of_distance_matrices <- map(list_of_distance_matrices, as.matrix)

# Training data is the source data
Dat.train <- LostInTheForest |> filter(class != "Human") |> droplevels() |> mutate(across(everything(), factor)) 

# Test data is the human data to be attributed to a source
Dat.test <- LostInTheForest |> filter(class == "Human") |> droplevels() |> mutate(across(everything(), factor)) 



## Source attribution of the sample data using the different methods
# prepare data using method of choice:
#' reminder of terminology: k = number ca axes (default is num.classes-1), 
#'                          m = number pco axes (default is num.classes-1), 
#'                          mp = propG (default is 100%), 
#'                          c = number cap axes (default is num.classes-1)
#'                          kc = number ca axes within cap (default is num.classes-1)

# 1. ca method
train <- prepare_training_ca(Dat.train, starts_with("Var"), "class")
test <- prepare_test_ca(Dat.test, train$extra, "id") 

# 2. binary method
train <- prepare_training_binary(Dat.train, starts_with("Var"), "class")
test <- prepare_test_binary(Dat.test, train$extra, "id")

# 3. ca unbiased method
train <- prepare_training_ca0(Dat.train, starts_with("Var"), "class")
test <- prepare_test_ca0(Dat.test, train$extra, "id")

# 4. pco method
train <- prepare_training_pco(Dat.train, starts_with("Var"), "class", d=list_of_distance_matrices, mp=99)
test <- prepare_test_pco(Dat.test, train$extra, "id")

# 5. cap method
train <- prepare_training_cap(Dat.train, starts_with("Var"), "class", d=list_of_distance_matrices, mp=99)
test <- prepare_test_cap(Dat.test, train$extra, "id")



## generate random forest models:
# 1. for binary method:
rf_mod <- ranger(dependent.variable.name = "class", data = train$training, oob.error = TRUE, num.trees=500, respect.unordered.factors = TRUE)

# 2. for other methods:
rf_mod <- ranger(class ~ ., data=train$training, oob.error = TRUE, num.trees=500, respect.unordered.factors = TRUE)



## make predictions for Human data
Prediction <- predict(rf_mod, data=test, predict.all = FALSE)$predictions
table(Prediction) |> as.data.frame() # counts of predictions for each source

## pull out individual tree decisions for each observation
uniques <- is_unique(Dat.test, train$extra)
tree_preds <- predict_by_tree(rf_mod, test, uniques)
answer <- tree_preds |> left_join(Dat.test |> select(id, class))
answer


if(0){
# test with variable with single level in training
Dat.train <- LostInTheForest %>% filter(class != "Human") |> mutate(Var01 = factor(10)) %>% droplevels() %>% mutate(across(everything(), factor)) 
train <- prepare_training_pco(Dat.train, starts_with("Var"), "class", d=list_of_distance_matrices, axes=1)
test <- prepare_test_pco(Dat.test, train$extra, "id")
}
