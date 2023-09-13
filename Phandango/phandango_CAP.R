### Run CAP with Phandango scores and compare to straight Hamming ###

## Load libraries and functions ##
source("methods/libs_fns.R")
source("methods/misclassification.R")
library(caret)   # for createFolds

## Load data ##
# allelic information - jejuni and coli 
load("../CAP_Data/data/cgMLST_dat.RData") # SACNZ cgMLST data set (cgMLST)

# distance information (Hamming)
load("../CAP_Data/data/list_of_distance_matrices.RData") # hamming distances

# distance information (Phandango)
load("Phandango/Phandango_data/list_of_distance_matrices_Phandango.RData") # hamming distances scaled to account for recombination

# species and clonal complex data
details <- read.csv("../CAP_data/data/SACNZ_referencelist.csv") |> select(LabID, Species, CC)

# prepare data - add details so can refine function later if required
Dat_jc <- cgMLST |> 
  left_join(details, by="LabID") |> select(LabID, Source, Species, CC, starts_with("CAMP")) |> 
  mutate(across(everything(), factor)) |> 
  mutate(Source = factor(Source, levels = c("Beef", "Poultry", "Sheep", "Human"))) |>
  mutate(Source = fct_recode(as.factor(Source), Cattle="Beef", Chicken="Poultry")) |> 
  filter(Source != "Human") |> droplevels()

# Create data splits
set.seed(3)
flds <- createFolds(y=Dat_jc$Source, k=10) 
Dat.train <- map(flds, ~slice(Dat_jc, {-.}))  # list of k training data sets
Dat.test <- map(flds, ~slice(Dat_jc, .))    # list of k test data sets


## 1. Test overall misclassification rates - cgMLST  # run 19/04/2023 - incorporated into "../CAP_Data/results/MC_results.Rdata"
MC_CAP.Ph <- misclass_fn(Dat.train, Dat.test, method="cap", d=list_of_distance_matrices_Phandango, mp=95, ntrees=500)
#load("../CAP_Data/results/MC_results.Rdata")

## 2. Test individual tree prediction accuracy # run 20/04/23 - incorporated into "../CAP_Data/results/results_cgMLST.Rdata"
cgMLST_CAP_P <- map2(Dat.train, Dat.test, ~tree_fn(.x,.y, method="cap", d=list_of_distance_matrices_Phandango, mp=95, class="Source", id="LabID", var_id="CAMP"))
#load("../CAP_Data/results/results_cgMLST.Rdata")

MC_all_trees <- misclass_tree_fn(results_cgMLST, class="Source")
#save(MC_all_trees, file="../CAP_Data/results/MC_all_trees.Rdata")  # updated 20/04/23

