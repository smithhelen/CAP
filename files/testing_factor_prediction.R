## Load libraries and functions ## 
source("methods/libs_fns.R") 
source("methods/misclassification.R") 
source("methods/tree_predictions.R") 

## Load data ## 
load("../CAP_data/data/cgMLST_dat.RData") # SACNZ cgMLST data set (jejuni and coli) 
details <- read.csv("../CAP_data/data/SACNZ_referencelist.csv") |> select(LabID, Species, CC) 
load("../CAP_data/data/CC_distance_matrices.Rdata")  # residualised based on CC 

# prepare data 
dat <- cgMLST |>  
  left_join(details, by="LabID") |> select(LabID, Source, Species, CC, starts_with("CAMP")) |>  
  mutate(across(everything(), factor)) |>  
  mutate(Source = factor(Source, levels = c("Beef", "Poultry", "Sheep", "Human"))) |> 
  mutate(Source = fct_recode(as.factor(Source), Cattle="Beef", Chicken="Poultry")) |>  
  filter(Source != "Human") |> droplevels() |> filter(!is.na(CC)) |> droplevels() |> arrange("LabID") # for residualised matrices based on clonal complex 
d <- map(CC_distance_matrices, as.matrix) # for residualised matrices based on clonal complex 
set.seed(3) 
samp <-sample(1:nrow(dat), 50, replace = FALSE) 
dat.train <- slice(dat, -samp) #|> select(1:6) 
dat.test <- slice(dat, samp) #|> select(1:6) 
alleles <- dat.train |> select(starts_with("CAMP")) |> colnames() 
d <- d[alleles]   

# residualisation
train <- prepare_training_pco(dat.train, starts_with("CAMP"), class="Source", d=d, axes=2, residualised="CC") 
test <- prepare_test_pco(dat.test, train$extra, id="LabID", residualised="CC") 

uniques <- is_unique(dat.test, train$extra)
set.seed(3) 
ranger_mod <- ranger(Source ~ ., data=train$training, oob.error = TRUE, num.trees=500, respect.unordered.factors = "partition", always.split.variables = "CC") 

preds_ranger <-predict(ranger_mod, data=test, predict.all = FALSE)$predictions 
(df <- data.frame(id = test |> pull("LabID"), predictions = preds_ranger) |> left_join(dat.test |> rename(id = "LabID") |> select(id,Source))) 

tree_preds_ranger <-predict(ranger_mod, data=test, predict.all = TRUE)$predictions 
(df_tree <- tree_preds_ranger |> as.data.frame() |>  
    mutate(id = test |> pull(LabID)) |> left_join(dat.test |> rename(id = LabID) |> select(id, Source)) |>  
    pivot_longer(-c(id, Source), names_to="tree", values_to="values") |>  
    left_join(data.frame(values = seq_along(ranger_mod$forest$levels), prediction = ranger_mod$forest$levels)) |>   
    mutate(tree = as.numeric(substring(tree, 2))) |>  
    select(id, tree, Source, prediction) |> arrange(tree)) 

tree_preds_JM <- predict_by_tree(ranger_mod, test, uniques, id="LabID", residualised="CC") 
(df_JM <- tree_preds_JM |> left_join(dat.test |> rename(id = LabID) |> select(id, Source)) |> select(id, tree, Source, prediction)) 

# check these are the same
df_JM |> anti_join(df_tree)

# YAY!

# no residualisation
train <- prepare_training_pco(dat.train, starts_with("CAMP"), class="Source", d=d, axes=2) 
test <- prepare_test_pco(dat.test, train$extra, id="LabID") 
uniques <- is_unique(dat.test, train$extra)
set.seed(3) 
ranger_mod <- ranger(Source ~ ., data=train$training, oob.error = TRUE, num.trees=500, respect.unordered.factors = TRUE)
preds_ranger <-predict(ranger_mod, data=test, predict.all = FALSE)$predictions 
(df <- data.frame(id = test |> pull("LabID"), predictions = preds_ranger) |> left_join(dat.test |> rename(id = "LabID") |> select(id,Source))) 
tree_preds_ranger <-predict(ranger_mod, data=test, predict.all = TRUE)$predictions 
(df_tree <- tree_preds_ranger |> as.data.frame() |>  
    mutate(id = test |> pull(LabID)) |> left_join(dat.test |> rename(id = LabID) |> select(id, Source)) |>  
    pivot_longer(-c(id, Source), names_to="tree", values_to="values") |>  
    left_join(data.frame(values = seq_along(ranger_mod$forest$levels), prediction = ranger_mod$forest$levels)) |>   
    mutate(tree = as.numeric(substring(tree, 2))) |>  
    select(id, tree, Source, prediction) |> arrange(tree)) 
tree_preds_JM <- predict_by_tree(ranger_mod, test, uniques, id="LabID", residualised=NULL) 
(df_JM <- tree_preds_JM |> left_join(dat.test |> rename(id = LabID) |> select(id, Source)) |> select(id, tree, Source, prediction)) 

# check these are the same
df_JM |> anti_join(df_tree)

# YAY!

