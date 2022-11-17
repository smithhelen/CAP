## Load libraries and functions ## 
source("methods/libs_fns.R") 
source("methods/misclassification.R") 
source("methods/tree_predictions.R") 

## Load data ## 
load("cgMLST_dat.RData") # SACNZ cgMLST data set (jejuni and coli) 
details <- read.csv("SACNZ_referencelist.csv") |> select(LabID, Species, CC) 
load("CC_distance_matrices.Rdata")  # residualised based on CC 

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


nd <-split(test, 1:nrow(test))  # list, each entry is a row of test data 
nu <-split(uniques, 1:nrow(uniques))  # list, each entry is a row of uniques (ie TRUE or FALSE for each var) 
id = test |> pull(LabID) 

fixup <- data.frame(prediction_update = ranger_mod$forest$levels[ranger_mod$forest$class.values],  
                    prediction_treeInfo = ranger_mod$forest$levels) 
tree1 <- treeInfo(ranger_mod,1) |> rename(prediction_treeInfo = prediction) |>  
  left_join(fixup, by="prediction_treeInfo") |> select(-prediction_treeInfo) |> rename(prediction = prediction_update) 

### tree1 
#  LabID   CC CAMP0001.V1 CAMP0001.V2 CAMP0002.V1 CAMP0002.V2 
# SC1439  828      1.4697   -1.640579   0.8025164   -1.143052 
predict_row(tree1, nd[[1]], nu[[1]], residualised="CC") # cattle 

#  LabID   CC CAMP0001.V1 CAMP0001.V2 CAMP0002.V1 CAMP0002.V2 
# SC0864  354    1.025315  -0.5974589   0.5343696    -0.95259 
predict_row(tree1, nd[[2]], nu[[2]], residualised="CC") # cattle 

#  LabID  CC CAMP0001.V1 CAMP0001.V2 CAMP0002.V1 CAMP0002.V2 
# SC1132  21    -2.58203   -8.437356    3.924993   -2.824844 
predict_row(tree1, nd[[3]], nu[[3]], residualised="CC") # chicken 

# with write.forest = TRUE 
ranger_mod2 <- ranger(Source ~ ., data=train$training, oob.error = TRUE, num.trees=5, respect.unordered.factors = "partition", always.split.variables = "CC", write.forest = TRUE) 

preds_ranger2 <-predict(ranger_mod2, data=test, predict.all = FALSE)$predictions 
(df2 <- data.frame(id = test |> pull("LabID"), predictions = preds_ranger2) |> left_join(dat.test |> rename(id = "LabID") |> select(id,Source))) 

tree_preds_ranger2 <-predict(ranger_mod2, data=test, predict.all = TRUE)$predictions 
(df_tree2 <- tree_preds_ranger2 |> as.data.frame() |>  
    mutate(id = test |> pull(LabID)) |> left_join(dat.test |> rename(id = LabID) |> select(id, Source)) |>  
    pivot_longer(-c(id, Source), names_to="tree", values_to="values") |>  
    left_join(data.frame(values = seq_along(ranger_mod2$forest$levels), prediction = ranger_mod2$forest$levels)) |>   
    mutate(tree = as.numeric(substring(tree, 2))) |>  
    select(id, tree, Source, prediction) |> arrange(tree)) 

tree_preds_JM2 <- predict_by_tree(ranger_mod2, test, uniques, id="LabID", residualised="CC") 
(df_JM2 <- tree_preds_JM2 |> left_join(dat.test |> rename(id = LabID) |> select(id, Source)) |> select(id, tree, Source, prediction)) 

fixup2 <- data.frame(prediction_update = ranger_mod2$forest$levels[ranger_mod2$forest$class.values],  
                     prediction_treeInfo = ranger_mod2$forest$levels) 
tree1b <- treeInfo(ranger_mod2,1) |> rename(prediction_treeInfo = prediction) |>  
  left_join(fixup2, by="prediction_treeInfo") |> select(-prediction_treeInfo) |> rename(prediction = prediction_update) 

all_equal(tree1, tree1b) 
all_equal(df, df2) 
all_equal(df_tree, df_tree2) 
all_equal(df_JM, df_JM2) 
