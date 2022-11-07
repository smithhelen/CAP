#### Split jejuni and coli
source("methods/libs_fns.R")
source("methods/tree_predictions.R")

## Load data - jejuni and coli combined
details <- read.csv("../CAP_data/data/SACNZ_referencelist.csv") |> select(LabID, Species)

load("../CAP_Data/data/cgMLST_dat.RData") # SACNZ cgMLST data set (jejuni and coli)
Dat_jc <- cgMLST |> mutate(across(everything(), factor)) |> 
  mutate(Source = factor(Source, levels = c("Beef", "Poultry", "Sheep", "Human"))) |>
  mutate(Source = fct_recode(as.factor(Source), Cattle="Beef", Chicken="Poultry")) |> 
  filter(Source != "Human") |> droplevels()

load("../CAP_Data/data/list_of_distance_matrices.RData") #  distance matrices of hamming distances among unique alleles for each gene
list_of_distance_matrices <- map(list_of_distance_matrices, as.matrix)

genes <- Dat_jc |> select(starts_with("CAMP")) |> colnames()

## function to subset distance matrix based on alleles ?not needed?
split_jc <- function(x,y, data){
  alleles <- data |> pull({{x}}) |> levels()
  d <- y[alleles,alleles]
  d
}


## Split data
# Training data is the jejuni data; test data is the coli data ###
# jejuni data
Dat_j <- Dat_jc |> left_join(details, by="LabID") |> filter(Species == "Jejuni") |> select(LabID, Source, Species, 2:1344) |> droplevels()
#dist_j <- map2(genes, list_of_distance_matrices, split_jc, Dat_j)

# coli data
Dat_c <- Dat_jc |> left_join(details, by="LabID") |> filter(Species == "Coli") |> select(LabID, Source, Species, 2:1344) |> droplevels()
#dist_c <- map2(genes, list_of_distance_matrices, split_jc, Dat_c)


## Run methods
#CA
ca.train <- prepare_training_ca0(Dat_j, starts_with("CAMP"), "Source", axes=2)
ca.test <- prepare_test_ca0(Dat_c, ca.train$extra, "LabID")
set.seed(3)
ca.rf_mod <- ranger(Source ~ ., data=ca.train$training, oob.error = TRUE, num.trees=500, respect.unordered.factors = TRUE)
predict(ca.rf_mod, data=ca.test, predict.all = FALSE)$predictions |> table(Dat_c$Source)
uniques.ca <- is_unique(Dat_c, ca.train$extra)
tree_preds.ca <- predict_by_tree(ca.rf_mod, ca.test, uniques.ca, "LabID")
answer.ca <- tree_preds.ca %>% left_join(Dat_c %>% rename(row = LabID) %>% select(row, Source))
answer.ca

#PCO
pco.train <- prepare_training_pco(Dat_j, starts_with("CAMP"), "Source", list_of_distance_matrices, axes=2)
pco.test <- prepare_test_pco(Dat_c, pco.train$extra, "LabID")
set.seed(3)
pco.rf_mod <- ranger(Source ~ ., data=pco.train$training, oob.error = TRUE, num.trees=500, respect.unordered.factors = TRUE)
predict(pco.rf_mod, data=pco.test, predict.all = FALSE)$predictions |> table(Dat_c$Source)
uniques.pco <- is_unique(Dat_c, pco.train$extra)
tree_preds.pco <- predict_by_tree(pco.rf_mod, pco.test, uniques.pco, "LabID")
answer.pco <- tree_preds.pco %>% left_join(Dat_c %>% rename(row = LabID) %>% select(row, Source))
answer.pco

#CAP
cap.train.99 <- prepare_training_cap(Dat_j, starts_with("CAMP"), "Source", list_of_distance_matrices, k=2, m=NULL, mp=99, axes=2)
cap.test.99 <- prepare_test_cap(Dat_c, cap.train.99$extra, "LabID")
set.seed(3)
cap.rf_mod.99 <- ranger(Source ~ ., data=cap.train.99$training, oob.error = TRUE, num.trees=500, respect.unordered.factors = TRUE)
predict(cap.rf_mod.99, data=cap.test.99, predict.all = FALSE)$predictions |> table(Dat_c$Source)
uniques.cap.99 <- is_unique(Dat_c, cap.train.99$extra)
tree_preds.cap.99 <- predict_by_tree(cap.rf_mod.99, cap.test.99, uniques.cap.99, "LabID")
answer.cap.99 <- tree_preds.cap.99 %>% left_join(Dat_c %>% rename(row = LabID) %>% select(row, Source))
answer.cap.99

cap.train.80 <- prepare_training_cap(Dat_j, starts_with("CAMP"), "Source", list_of_distance_matrices, k=2, m=NULL, mp=80, axes=2)
cap.test.80 <- prepare_test_cap(Dat_c, cap.train.80$extra, "LabID")
set.seed(3)
cap.rf_mod.80 <- ranger(Source ~ ., data=cap.train.80$training, oob.error = TRUE, num.trees=500, respect.unordered.factors = TRUE)
predict(cap.rf_mod.80, data=cap.test.80, predict.all = FALSE)$predictions |> table(Dat_c$Source)
uniques.cap.80 <- is_unique(Dat_c, cap.train.80$extra)
tree_preds.cap.80 <- predict_by_tree(cap.rf_mod.80, cap.test.80, uniques.cap.80, "LabID")
answer.cap.80 <- tree_preds.cap.80 %>% left_join(Dat_c %>% rename(row = LabID) %>% select(row, Source))
answer.cap.80

