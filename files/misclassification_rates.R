# pull out individual tree predictions and splitting genes for cgMLST and identify when unique alleles were used
# these take a long time as they run through the 10 folds of 500 trees
# compare for CAP_95 and CAP_CC_95 only

## Load libraries and functions #############
source("methods/libs_fns.R")
source("methods/tree_predictions.R")  # to pull out individual tree predictions
source("methods/misclassification.R")  # to calculate mc rates
source("methods/ranger_mods.R")
library(caret)   # for createFolds


# define function to prepare training and test data, train the rf model, and calculate misclassification rates (of trees and forests)
get_mc <- function(train.flds, test.flds, method="ca0", residualised=NULL, d=NULL, k=NULL, m=NULL, mp=99, c=NULL, id, class, var_id, ntrees=500, tree=TRUE, forest=TRUE){
  model.info.flds <- map2(train.flds, test.flds, function(Dat.train, Dat.test) {
    prep_and_train(Dat.train, Dat.test, method=method, residualised, d=d, k=k, m=m, mp=mp, c=c, id=id, class=class, var_id=var_id, ntrees=ntrees)})
  results <- mc_function(test.flds, model.info.flds, class=class, tree=tree, forest=forest)
  results
}



## Load data #############
# allelic information - jejuni and coli 
load("../CAP_Data/data/cgMLST_dat.RData") # SACNZ cgMLST data set (jejuni and coli)

# species and clonal complex data
details <- read.csv("../CAP_data/data/SACNZ_referencelist.csv") |> select(LabID, Species, CC)

# distance information (Hamming)
load("../CAP_Data/data/list_of_distance_matrices.RData") #  distance matrices of hamming distances among unique alleles for each gene

# distance information (Phandango)
load("Phandango/Phandango_data/list_of_distance_matrices_Phandango.RData") # hamming distances scaled to account for recombination

# residualised distance information
load("../CAP_Data/data/CC_distance_matrices.Rdata")  # residualised based on CC

## Extra data ##
# amino acid alleleic information
load("../CAP_data/data/aa_allele_dat.RData")

# residualised distance information
load("../CAP_Data/data/species_distance_matrices.Rdata") #  residualised based on species

# amino acid distance information
load("../CAP_data/data/aa_distance_matrices.RData")
load("../CAP_Data/data/CC_distance_matrices_aa.Rdata")



## Select dataset to run cross-validation on #############
Dat <- cgMLST |> 
  left_join(details, by="LabID") |> select(LabID, Source, Species, CC, starts_with("CAMP")) |> 
  mutate(across(everything(), factor)) |> 
  mutate(Source = factor(Source, levels = c("Beef", "Poultry", "Sheep", "Human"))) |>
  mutate(Source = fct_recode(as.factor(Source), Cattle="Beef", Chicken="Poultry")) |> 
  filter(Source != "Human") |> droplevels()
Dat_aa <- aa_allele_dat |> 
  left_join(details, by="LabID") |> select(LabID, Source, Species, CC, starts_with("CAMP")) |> 
  mutate(across(everything(), factor)) |> 
  mutate(Source = factor(Source, levels = c("Beef", "Poultry", "Sheep", "Human"))) |>
  mutate(Source = fct_recode(as.factor(Source), Cattle="Beef", Chicken="Poultry")) |> 
  filter(Source != "Human") |> droplevels()

dat <- Dat # original jejuni coli data
dat_CC <- Dat |> filter(!is.na(CC)) |> droplevels() |> arrange("LabID") # for residualised matrices based on clonal complex
dat_sp <- Dat |> filter(!is.na(Species)) |> droplevels() # for residualised matrices based on species
dat_j <- Dat |> filter(Species=="Jejuni") |> droplevels() # jejuni only data
dat_aa <- Dat_aa
dat_aa_CC <- Dat_aa|> filter(!is.na(CC)) |> droplevels() |> arrange("LabID") # for aa residualised data

## Split data into training and test sets using folds #############
# 1. original data
set.seed(3)
flds <- createFolds(y=dat$Source, k=10) 
Dat.train <- map(flds, ~slice(dat, {-.}))  # list of k training data sets
Dat.test <- map(flds, ~slice(dat, .))    # list of k test data sets

# 2. CC residualised data
set.seed(3)
flds.CC <- createFolds(y=dat_CC$Source, k=10) 
Dat.train.CC <- map(flds.CC, ~slice(dat_CC, {-.}))  # list of k training data sets
Dat.test.CC <- map(flds.CC, ~slice(dat_CC, .))    # list of k test data sets

# 3. species residualised data
set.seed(3)
flds.sp <- createFolds(y=dat_sp$Source, k=10) 
Dat.train.sp <- map(flds.sp, ~slice(dat_sp, {-.}))  # list of k training data sets
Dat.test.sp <- map(flds.sp, ~slice(dat_sp, .))    # list of k test data sets

# 4. jejuni only data
set.seed(3)
flds.j <- createFolds(y=dat_j$Source, k=10) 
Dat.train.j <- map(flds.j, ~slice(dat_j, {-.}))  # list of k training data sets
Dat.test.j <- map(flds.j, ~slice(dat_j, .))    # list of k test data sets

# 5. amino acid data
set.seed(3)
flds.aa <- createFolds(y=dat_aa$Source, k=10) 
Dat.train.aa <- map(flds.aa, ~slice(dat_aa, {-.}))  # list of k training data sets
Dat.test.aa <- map(flds.aa, ~slice(dat_aa, .))    # list of k test data sets

# 6. amino acid CC residualised data
set.seed(3)
flds.aa_CC <- createFolds(y=dat_aa_CC$Source, k=10) 
Dat.train.aa_CC <- map(flds.aa_CC, ~slice(dat_aa_CC, {-.}))  # list of k training data sets
Dat.test.aa_CC <- map(flds.aa_CC, ~slice(dat_aa_CC, .))    # list of k test data sets



## Select distance matrices to match choice of data #################
d <- map(list_of_distance_matrices, as.matrix) # original
d_CC <- map(CC_distance_matrices, as.matrix) # for residualised matrices based on clonal complex
d_P <- map(list_of_distance_matrices_Phandango, as.matrix) # Phandango (recombination)
d_sp <- map(species_distance_matrices, as.matrix) # for residualised matrices based on species
d_aa <- map(aa_distance_matrices, as.matrix) # amino acid data
d_aa_CC <- map(CC_distance_matrices_aa, as.matrix) # for amino acid residualised matrices based on clonal complex



## Pull out forest and tree predictions #############
doParallel::registerDoParallel()

#' reminder of terminology: k = number ca axes (default is num.classes-1), 
#'                          m = number pco axes (default is num.classes-1), 
#'                          mp = propG (default is 100%), 
#'                          c = number cap axes (default is num.classes-1)

# 1. original data, 1 axis, standard distance measure
cgMLST_CA0_1 <- get_mc(Dat.train, Dat.test, method="ca0", class="Source", id="LabID", var_id= "CAMP", k=1)
cgMLST_PCO_1 <- get_mc(Dat.train, Dat.test, method="pco", class="Source", id="LabID", var_id="CAMP", d=d, m=1)
cgMLST_CAP_1 <- get_mc(Dat.train, Dat.test, method="cap", class="Source", id="LabID", var_id="CAMP", d=d, k=2, mp=95, c=1)

# 1b. original data, mp=95 or 2 axes, standard distance measure
cgMLST_CA0 <- get_mc(Dat.train, Dat.test, method="ca0", class="Source", id="LabID", var_id="CAMP", k=2)
cgMLST_PCO <- get_mc(Dat.train, Dat.test, method="pco", class="Source", id="LabID", var_id="CAMP", d=d, mp=95)
cgMLST_CAP <- get_mc(Dat.train, Dat.test, method="cap", class="Source", id="LabID", var_id="CAMP", d=d, k=2, mp=95, c=2)
cgMLST_CA0_2 <- cgMLST_CA0
cgMLST_PCO_2 <- cgMLST_PCO
cgMLST_CAP_2 <- cgMLST_CAP

# 1c. original data, mp=95, Phandango distance measure
cgMLST_PCO_Ph <- get_mc(Dat.train, Dat.test, method="pco", class="Source", id="LabID", var_id="CAMP", d=d_P, mp=95)
cgMLST_CAP_Ph <- get_mc(Dat.train, Dat.test, method="cap", class="Source", id="LabID", var_id="CAMP", d=d_P, k=2, mp=95, c=2)

# 1c. original data, 2 cap axes, change mp values
cgMLST_CAP_0 <- get_mc(Dat.train, Dat.test, method="cap", class="Source", id="LabID", var_id="CAMP", d=d, k=2, mp=0, c=2)
cgMLST_CAP_50 <- get_mc(Dat.train, Dat.test, method="cap", class="Source", id="LabID", var_id="CAMP", d=d, k=2, mp=50, c=2)
cgMLST_CAP_85 <- get_mc(Dat.train, Dat.test, method="cap", class="Source", id="LabID", var_id="CAMP", d=d, k=2, mp=85, c=2)
cgMLST_CAP_95 <- cgMLST_CAP
cgMLST_CAP_99 <- get_mc(Dat.train, Dat.test, method="cap", class="Source", id="LabID", var_id="CAMP", d=d, k=2, mp=99, c=2)

# 2. CC residualised data, mp=95 or 2 axes, CC distance measure
cgMLST_CA0_CC <- get_mc(Dat.train.CC, Dat.test.CC, method="ca0", class="Source", id="LabID", var_id="CAMP", residualised="CC", k=2)
cgMLST_PCO_CC <- get_mc(Dat.train.CC, Dat.test.CC, method="pco", class="Source", id="LabID", var_id="CAMP", residualised="CC", d=d_CC, mp=95)
cgMLST_CAP_CC <- get_mc(Dat.train.CC, Dat.test.CC, method="cap", class="Source", id="LabID", var_id="CAMP", residualised="CC", d=d_CC, k=2, mp=95, c=2)

# 3. species residualised data, mp=95 or 2 axes, species distance measure
cgMLST_CA0_sp <- get_mc(Dat.train.sp, Dat.test.sp, method="ca0", class="Source", id="LabID", var_id="CAMP", residualised="Species", k=2)
cgMLST_PCO_sp <- get_mc(Dat.train.sp, Dat.test.sp, method="pco", class="Source", id="LabID", var_id="CAMP", residualised="Species", d=d_sp, mp=95)
cgMLST_CAP_sp <- get_mc(Dat.train.sp, Dat.test.sp, method="cap", class="Source", id="LabID", var_id="CAMP", residualised="Species", d=d_sp, k=2, mp=95, c=2)

# 4. jejuni only data, mp=95 or 2 axes, original distance measure
cgMLST_CA0_j <- get_mc(Dat.train.j, Dat.test.j, method="ca0", class="Source", id="LabID", var_id="CAMP", k=2)
cgMLST_PCO_j <- get_mc(Dat.train.j, Dat.test.j, method="pco", class="Source", id="LabID", var_id="CAMP", d=d, mp=95)
cgMLST_CAP_j <- get_mc(Dat.train.j, Dat.test.j, method="cap", class="Source", id="LabID", var_id="CAMP", d=d, k=2, mp=95, c=2)

# 5. amino acid data, mp=95 or 2 axes, aa distance measure
cgMLST_CA0_aa <- get_mc(Dat.train.aa, Dat.test.aa, method="ca0", class="Source", id="LabID", var_id="CAMP", k=2)
cgMLST_PCO_aa <- get_mc(Dat.train.aa, Dat.test.aa, method="pco", class="Source", id="LabID", var_id="CAMP", d=d_aa, mp=95)
cgMLST_CAP_aa <- get_mc(Dat.train.aa, Dat.test.aa, method="cap", class="Source", id="LabID", var_id="CAMP", d=d_aa, k=2, mp=95, c=2)

# 6. amino acid CC residualised data, mp=95 or 2 axes, aa_CC distance measure
cgMLST_CA0_aa_CC <- get_mc(Dat.train.aa_CC, Dat.test.aa_CC, method="ca0", class="Source", id="LabID", var_id="CAMP", residualised="CC", k=2)
cgMLST_PCO_aa_CC <- get_mc(Dat.train.aa_CC, Dat.test.aa_CC, method="pco", class="Source", id="LabID", var_id="CAMP", residualised="CC", d=d_aa_CC, mp=95)
cgMLST_CAP_aa_CC <- get_mc(Dat.train.aa_CC, Dat.test.aa_CC, method="cap", class="Source", id="LabID", var_id="CAMP", residualised="CC", d=d_aa_CC, k=2, mp=95, c=2)


## Merge data sets
cgMLST_mc_results <- bind_rows(CA0=cgMLST_CA0, PCO=cgMLST_PCO, CAP=cgMLST_CAP, 
                               CA0_CC=cgMLST_CA0_CC, PCO_CC=cgMLST_PCO_CC, CAP_CC=cgMLST_CAP_CC,
                               PCO_Ph=cgMLST_CAP_Ph, CAP_Ph=cgMLST_CAP_Ph,
                               CA01=cgMLST_CA0_1, PCO1=cgMLST_PCO_1, CAP1=cgMLST_CAP_1,
                               CA02=cgMLST_CA0, PCO2=cgMLST_PCO, CAP2=cgMLST_CAP,
                               CAP_0=cgMLST_CAP_0, CAP_50=cgMLST_CAP_50, CAP_80=cgMLST_CAP_85, CAP_95=cgMLST_CAP, CAP_99=cgMLST_CAP_99,
                               CA0_sp=cgMLST_CA0_sp, PCO_sp=cgMLST_PCO_sp, CAP_sp=cgMLST_CAP_sp,
                               CA0_j=cgMLST_CA0_j, PCO_j=cgMLST_PCO_j, CAP_j=cgMLST_CAP_j,
                               CA0_aa=cgMLST_CA0_aa, PCO_aa=cgMLST_PCO_aa, CAP_aa=cgMLST_CAP_aa,
                               CA0_aa_CC=cgMLST_CA0_aa_CC, PCO_aa_CC=cgMLST_PCO_aa_CC, CAP_aa_CC=cgMLST_CAP_aa_CC,
                               .id="method")
#save(cgMLST_mc_results, file = "../CAP_Data/results/cgMLST_mc_results.Rdata")  # this was run on 16 June


# now can just load :-)
load("../CAP_Data/results/cgMLST_mc_results.Rdata")

# and plot


# to get fold information, need to re-run but select tree=FALSE







