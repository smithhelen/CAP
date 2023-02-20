# pull out individual tree predictions and splitting genes for cgMLST and identify when unique alleles were used
# these take a long time as they run through the 10 folds of 500 trees
# compare for CAP_95 and CAP_CC_95 only

## Load libraries and functions
source("methods/libs_fns.R")
source("methods/tree_predictions.R")  # to pull out individual tree predictions
library(caret)   # for createFolds


## Load data
# allelic information - jejuni and coli 
load("../CAP_Data/data/cgMLST_dat.RData") # SACNZ cgMLST data set (jejuni and coli)

# species and clonal complex data
details <- read.csv("../CAP_data/data/SACNZ_referencelist.csv") |> select(LabID, Species, CC)

# distance information (Hamming)
load("../CAP_Data/data/list_of_distance_matrices.RData") #  distance matrices of hamming distances among unique alleles for each gene

# residualised distance information
load("../CAP_Data/data/CC_distance_matrices.Rdata")  # residualised based on CC


## Prepare data
Dat <- cgMLST |> 
  left_join(details, by="LabID") |> select(LabID, Source, Species, CC, starts_with("CAMP")) |> 
  mutate(across(everything(), factor)) |> 
  mutate(Source = factor(Source, levels = c("Beef", "Poultry", "Sheep", "Human"))) |>
  mutate(Source = fct_recode(as.factor(Source), Cattle="Beef", Chicken="Poultry")) |> 
  filter(Source != "Human") |> droplevels()

# Select dataset to run cross-validation on
dat_95 <- Dat # original jejuni coli data
dat_CC <- Dat |> filter(!is.na(CC)) |> droplevels() |> arrange("LabID") # for residualised matrices based on clonal complex

# Select distance matrices to match choice of data
d_95 <- map(list_of_distance_matrices, as.matrix) # original
d_CC <- map(CC_distance_matrices, as.matrix) # for residualised matrices based on clonal complex


## Create data splits - repeat for each different set of data being run ##
# 1. original data, mp=95
set.seed(3)
flds <- createFolds(y=dat_95$Source, k=10) 
Dat.train <- map(flds, ~slice(dat_95, {-.}))  # list of k training data sets
Dat.test <- map(flds, ~slice(dat_95, .))    # list of k test data sets

# 2. CC residualised data, mp=95
set.seed(3)
flds.CC <- createFolds(y=dat_CC$Source, k=10) 
Dat.train.CC <- map(flds.CC, ~slice(dat_CC, {-.}))  # list of k training data sets
Dat.test.CC <- map(flds.CC, ~slice(dat_CC, .))    # list of k test data sets


## Pull out individual tree predictions
# 1. original data, mp=95
cgMLST_CA0 <- map2(Dat.train, Dat.test, ~tree_fn(.x,.y, method="ca0", axes=2, class="Source", id="LabID", var_id="CAMP"))
cgMLST_PCO <- map2(Dat.train, Dat.test, ~tree_fn(.x,.y, d=d_95, method="pco", axes=2, class="Source", id="LabID", var_id="CAMP"))
cgMLST_CAP <- map2(Dat.train, Dat.test, ~tree_fn(.x,.y, d=d_95, method="cap", k=2, mp=95, axes=1, class="Source", id="LabID", var_id="CAMP"))

# 2. CC residualised data, mp=95
cgMLST_CA0_CC <- map2(Dat.train.CC, Dat.test.CC, ~tree_fn(.x,.y, method="ca0", axes=2, class="Source", id="LabID", residualised="CC", var_id="CAMP"))
cgMLST_PCO_CC <- map2(Dat.train.CC, Dat.test.CC, ~tree_fn(.x,.y, d=d_CC, method="pco", class="Source", id="LabID", axes=2, residualised="CC", var_id="CAMP"))
cgMLST_CAP_CC <- map2(Dat.train.CC, Dat.test.CC, ~tree_fn(.x,.y, d=d_CC, method="cap", class="Source", id="LabID", k=2, mp=95, axes=2, residualised="CC", var_id="CAMP"))


# Merge data sets
results_cgMLST <- bind_rows(CA0=cgMLST_CA0 |> bind_rows(.id = "Fold"), 
                            PCO=cgMLST_PCO |> bind_rows(.id = "Fold"), 
                            CAP=cgMLST_CAP |> bind_rows(.id = "Fold"), 
                            CA0_CC=cgMLST_CA0_CC |> bind_rows(.id = "Fold"), 
                            PCO_CC=cgMLST_PCO_CC |> bind_rows(.id = "Fold"), 
                            CAP_CC=cgMLST_CAP_CC |> bind_rows(.id = "Fold"), 
                            .id="method")
save(results_cgMLST, file = "../CAP_Data/results/results_cgMLST.Rdata")  # this was run on 18/11/2022


# now can just load :-)
load("../CAP_Data/results/results_cgMLST.Rdata")

# and use for cross_validation.R



