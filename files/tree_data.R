# pull out individual tree predictions and splitting genes 
# for cgMLST and for MLST
# and identify when unique alleles were used
# these take a long time as they run through the 10 folds of 500 trees

# load libraries and functions
source("methods/libs_fns.R")
source("methods/tree_predictions.R")  # to pull out individual tree predictions

# load data - jejuni and coli combined
load("../CAP_Data/data/cgMLST_dat.RData") # SACNZ cgMLST data set (jejuni and coli)
Dat_jc <- cgMLST |> mutate(across(everything(), factor)) |> 
  mutate(Source = factor(Source, levels = c("Beef", "Poultry", "Sheep", "Human"))) |>
  mutate(Source = fct_recode(as.factor(Source), Cattle="Beef", Chicken="Poultry")) |> 
  filter(Source != "Human") |> droplevels()

# load distance matrices
load("../CAP_Data/data/list_of_distance_matrices.RData") #  distance matrices of hamming distances among unique alleles for each gene
list_of_distance_matrices <- map(list_of_distance_matrices, as.matrix)

# create data splits
set.seed(3)
flds <- createFolds(y=Dat_jc$Source, k=10) 
Dat.train <- map(flds, ~slice(Dat_jc, {-.}))  # list of k training data sets
Dat.test <- map(flds, ~slice(Dat_jc, .))    # list of k test data sets

## 1. Pull out individual tree predictions

cgMLST_CA.zero1_jc <- map2(Dat.train, Dat.test, ~tree_fn(.x,.y, method="ca0", axes=1))
cgMLST_CA.zero2_jc <- map2(Dat.train, Dat.test, ~tree_fn(.x,.y, method="ca0", axes=2))
cgMLST_PCO1_jc <- map2(Dat.train, Dat.test, ~tree_fn(.x,.y, d=list_of_distance_matrices, method="pco", axes=1))
cgMLST_PCO2_jc <- map2(Dat.train, Dat.test, ~tree_fn(.x,.y, d=list_of_distance_matrices, method="pco", axes=2))
cgMLST_CAP1_jc <- map2(Dat.train, Dat.test, ~tree_fn(.x,.y, d=list_of_distance_matrices, method="cap", k=2, m=2, axes=1))
cgMLST_CAP2_jc <- map2(Dat.train, Dat.test, ~tree_fn(.x,.y, d=list_of_distance_matrices, method="cap", k=2, m=2, axes=2))

# Merge data sets
results_cgMLST_jc <- bind_rows(CA1=cgMLST_CA.zero1_jc |> bind_rows(.id = "Fold"), 
                               CA2=cgMLST_CA.zero2_jc |> bind_rows(.id = "Fold"), 
                               PCO1=cgMLST_PCO1_jc |> bind_rows(.id = "Fold"), 
                               PCO2=cgMLST_PCO2_jc |> bind_rows(.id = "Fold"), 
                               CAP1=cgMLST_CAP1_jc |> bind_rows(.id = "Fold"), 
                               CAP2=cgMLST_CAP2_jc |> bind_rows(.id = "Fold"), 
                               .id="method")
save(results_cgMLST_jc, file = "../CAP_Data/results/results_cgMLST_jc.Rdata")  # this was run on ...

# now can just load :-)
load("../CAP_Data/results/results_cgMLST_jc.Rdata")

## 2. try different mp values to determine number of axes
#mp_list <- list(0,50,85,95,99)
mp_list <- list(95)
cgMLST_CAP_jc <- map(mp_list, function(z) map2(Dat.train, Dat.test, ~tree_fn(.x,.y,mp=z,d=list_of_distance_matrices, method="cap", k=2, axes=2)))

# check which answer_CAP_jc in the list is for mp=95
map(cgMLST_CAP_jc, ~.x |> map("mp") |> unlist() |> unique())

results_cgMLST_jc <- bind_rows(results_cgMLST_jc, CAP=cgMLST_CAP_jc[[1]] |> map("answer") |> bind_rows(.id = "Fold"), .id="method")
save(results_cgMLST_jc, file = "../CAP_Data/results/results_cgMLST_jc.Rdata")




