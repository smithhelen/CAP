## Useful for testing the code is running
## Source attribution of the sample data using the different methods

# Load libraries and functions
source("methods/libs_fns.R")

# Load data
load("../CAP_Data/data/cgMLST_dat.RData") # SACNZ cgMLST data set (jejuni and coli)
Dat_jc <- cgMLST |> mutate(across(everything(), factor)) |> 
  mutate(Source = factor(Source, levels = c("Beef", "Poultry", "Sheep", "Human"))) |>
  mutate(Source = fct_recode(as.factor(Source), Cattle="Beef", Chicken="Poultry")) 

load("../CAP_Data/data/list_of_distance_matrices_all.RData") #  distance matrices of hamming distances among unique alleles for each gene
list_of_distance_matrices <- map(list_of_distance_matrices_all, as.matrix)

# Subset data
small_dat <- cgMLST |> select(LabID, Source, sample(2:1344,10))
genes <- small_dat |> select(starts_with("CAMP")) |> colnames()
small_d <- list_of_distance_matrices[genes]

# Create data splits
Dat.train <- small_dat |> filter(Source != "Human") |> droplevels() 
Dat.test <- small_dat |> filter(Source == "Human") |> droplevels() 

# prepare data using method of choice:
# 1. ca unbiased method, 1 axis
train <- prepare_training_ca0(Dat.train, starts_with("CAMP"), "Source", axes=1)
test <- prepare_test_ca0(Dat.test, train$extra, "LabID")
rf_mod <- ranger(Source ~ ., data=train$training, oob.error = TRUE, num.trees=500, respect.unordered.factors = TRUE)
Prediction <- predict(rf_mod, data=test, predict.all = FALSE)$predictions
table(Prediction) |> as.data.frame() # counts of predictions for each source
#uniques <- is_unique(Dat.test, train$extra)
#tree_preds <- predict_by_tree(rf_mod, test, uniques)
#answer <- tree_preds |> left_join(Dat.test |> rename(row = id) |> select(row, class))
#answer

# 2. ca unbiased method, 2 axes
train <- prepare_training_ca0(Dat.train, starts_with("CAMP"), "Source", axes=2)
test <- prepare_test_ca0(Dat.test, train$extra, "LabID")
rf_mod <- ranger(Source ~ ., data=train$training, oob.error = TRUE, num.trees=500, respect.unordered.factors = TRUE)
Prediction <- predict(rf_mod, data=test, predict.all = FALSE)$predictions
table(Prediction) |> as.data.frame() # counts of predictions for each source
#uniques <- is_unique(Dat.test, train$extra)
#tree_preds <- predict_by_tree(rf_mod, test, uniques)
#answer <- tree_preds |> left_join(Dat.test |> rename(row = id) |> select(row, class))
#answer

# 3. pco method, 1 axis
train <- prepare_training_pco(Dat.train, starts_with("CAMP"), "Source", d=small_d, axes=1)
test <- prepare_test_pco(Dat.test, train$extra, "LabID")
rf_mod <- ranger(Source ~ ., data=train$training, oob.error = TRUE, num.trees=500, respect.unordered.factors = TRUE)
Prediction <- predict(rf_mod, data=test, predict.all = FALSE)$predictions
table(Prediction) |> as.data.frame() # counts of predictions for each source
#uniques <- is_unique(Dat.test, train$extra)
#tree_preds <- predict_by_tree(rf_mod, test, uniques)
#answer <- tree_preds |> left_join(Dat.test |> rename(row = id) |> select(row, class))
#answer

# 4. pco method, 2 axes
train <- prepare_training_pco(Dat.train, starts_with("CAMP"), "Source", d=small_d, axes=2)
test <- prepare_test_pco(Dat.test, train$extra, "LabID")
rf_mod <- ranger(Source ~ ., data=train$training, oob.error = TRUE, num.trees=500, respect.unordered.factors = TRUE)
Prediction <- predict(rf_mod, data=test, predict.all = FALSE)$predictions
table(Prediction) |> as.data.frame() # counts of predictions for each source
#uniques <- is_unique(Dat.test, train$extra)
#tree_preds <- predict_by_tree(rf_mod, test, uniques)
#answer <- tree_preds |> left_join(Dat.test |> rename(row = id) |> select(row, class))
#answer

# 5. cap method, 1 axis
train <- prepare_training_cap(Dat.train, starts_with("CAMP"), "Source", d=small_d, mp=95, axes=1)
test <- prepare_test_cap(Dat.test, train$extra, "LabID")
rf_mod <- ranger(Source ~ ., data=train$training, oob.error = TRUE, num.trees=500, respect.unordered.factors = TRUE)
Prediction <- predict(rf_mod, data=test, predict.all = FALSE)$predictions
table(Prediction) |> as.data.frame() # counts of predictions for each source
#uniques <- is_unique(Dat.test, train$extra)
#tree_preds <- predict_by_tree(rf_mod, test, uniques)
#answer <- tree_preds |> left_join(Dat.test |> rename(row = id) |> select(row, class))
#answer

# 6. cap method, 2 axes
train <- prepare_training_cap(Dat.train, starts_with("CAMP"), "Source", d=small_d, mp=95, axes=2)
test <- prepare_test_cap(Dat.test, train$extra, "LabID")
rf_mod <- ranger(Source ~ ., data=train$training, oob.error = TRUE, num.trees=500, respect.unordered.factors = TRUE)
Prediction <- predict(rf_mod, data=test, predict.all = FALSE)$predictions
table(Prediction) |> as.data.frame() # counts of predictions for each source
#uniques <- is_unique(Dat.test, train$extra)
#tree_preds <- predict_by_tree(rf_mod, test, uniques)
#answer <- tree_preds |> left_join(Dat.test |> rename(row = id) |> select(row, class))
#answer



