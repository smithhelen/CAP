## calculate (1) misclassification rates of the individual trees according to the presence of absent levels 
# and (2) as overal misclassification rates of the whole forests for the different methods
# Compare cgMLST and MLST

## Load libraries and functions ##
source("methods/libs_fns.R")
source("methods/misclassification.R")
library(caret)   # for createFolds

## Load data ##
# allelic information - jejuni and coli 
load("../CAP_Data/data/cgMLST_dat.RData") # SACNZ cgMLST data set (jejuni and coli)

# amino acid alleleic information
load("../CAP_data/data/aa_allele_dat.RData")

# species and clonal complex data
details <- read.csv("../CAP_data/data/SACNZ_referencelist.csv") |> select(LabID, Species, CC)

# distance information (Hamming)
load("../CAP_Data/data/list_of_distance_matrices.RData") #  distance matrices of hamming distances among unique alleles for each gene

# distance information (Phandango)
load("Phandango/Phandango_data/list_of_distance_matrices_Phandango.RData") # hamming distances scaled to account for recombination

# residualised distance information
load("../CAP_Data/data/species_distance_matrices.Rdata") #  residualised based on species
load("../CAP_Data/data/CC_distance_matrices.Rdata")  # residualised based on CC

# amino acid distance information
load("../CAP_data/data/aa_distance_matrices.RData")
load("../CAP_Data/data/CC_distance_matrices_aa.Rdata")

# prepare data
Dat_jc <- cgMLST |> 
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


# Select dataset to run cross-validation on
dat <- Dat_jc # original jejuni coli data
dat_sp <- Dat_jc |> filter(!is.na(Species)) |> droplevels() # for residualised matrices based on species
dat_CC <- Dat_jc |> filter(!is.na(CC)) |> droplevels() |> arrange("LabID") # for residualised matrices based on clonal complex
dat_j <- Dat_jc |> filter(Species=="Jejuni") |> droplevels() # jejuni only data
dat_aa <- Dat_aa
dat_aa_CC <- Dat_aa|> filter(!is.na(CC)) |> droplevels() |> arrange("LabID") # for aa residualised data

# Select distance matrices to match choice of data
list_of_distance_matrices <- map(list_of_distance_matrices, as.matrix) # original
list_of_distance_matrices_sp <- map(species_distance_matrices, as.matrix) # for residualised matrices based on species
list_of_distance_matrices_CC <- map(CC_distance_matrices, as.matrix) # for residualised matrices based on clonal complex
list_of_distance_matrices_aa <- map(aa_distance_matrices, as.matrix) # amino acid data
list_of_distance_matrices_aa_CC <- map(CC_distance_matrices_aa, as.matrix) # for amino acid residualised matrices based on clonal complex
list_of_distance_matrices_Phandango <- map(list_of_distance_matrices_Phandango, as.matrix) # Phandango (recombination)


## Create data splits - repeat for each different set of data being run ##
set.seed(3)
flds <- createFolds(y=dat$Source, k=10) 
Dat.train <- map(flds, ~slice(dat, {-.}))  # list of k training data sets
Dat.test <- map(flds, ~slice(dat, .))    # list of k test data sets

set.seed(3)
flds_sp <- createFolds(y=dat_sp$Source, k=10) 
Dat.train_sp <- map(flds_sp, ~slice(dat_sp, {-.}))  # list of k training data sets
Dat.test_sp <- map(flds_sp, ~slice(dat_sp, .))    # list of k test data sets

set.seed(3)
flds_CC <- createFolds(y=dat_CC$Source, k=10) 
Dat.train_CC <- map(flds_CC, ~slice(dat_CC, {-.}))  # list of k training data sets
Dat.test_CC <- map(flds_CC, ~slice(dat_CC, .))    # list of k test data sets

set.seed(3)
flds_j <- createFolds(y=dat_j$Source, k=10) 
Dat.train_j <- map(flds_j, ~slice(dat_j, {-.}))  # list of k training data sets
Dat.test_j <- map(flds_j, ~slice(dat_j, .))    # list of k test data sets

set.seed(3)
flds_aa <- createFolds(y=dat_aa$Source, k=10) 
Dat.train_aa <- map(flds_aa, ~slice(dat_aa, {-.}))  # list of k training data sets
Dat.test_aa <- map(flds_aa, ~slice(dat_aa, .))    # list of k test data sets

set.seed(3)
flds_aa_CC <- createFolds(y=dat_aa_CC$Source, k=10) 
Dat.train_aa_CC <- map(flds_aa_CC, ~slice(dat_aa_CC, {-.}))  # list of k training data sets
Dat.test_aa_CC <- map(flds_aa_CC, ~slice(dat_aa_CC, .))    # list of k test data sets


## 1. Test overall misclassification rates - cgMLST  # run 27/10/2022, fully re-run 20/04/23
doParallel::registerDoParallel()

# 1 axis 
MC_CA01 <- misclass_fn(Dat.train, Dat.test, method="ca0", axes=1)
MC_PCO1 <- misclass_fn(Dat.train, Dat.test, d=list_of_distance_matrices, method="pco", axes=1)
MC_CAP1 <- misclass_fn(Dat.train, Dat.test, d=list_of_distance_matrices, method="cap", k=2, mp=95, axes=1) # mp (prop. of variance explained by m PCO axes, aka propG)

# 2 axes
MC_CA02 <- misclass_fn(Dat.train, Dat.test, method="ca0", axes=2)
MC_PCO2 <- misclass_fn(Dat.train, Dat.test, d=list_of_distance_matrices, method="pco", axes=2)
MC_CAP2 <- misclass_fn(Dat.train, Dat.test, d=list_of_distance_matrices, method="cap", k=2, mp=95, axes=2)

# residualised matrices
# what about respect.unordered.factors="partition" (the misclassification code was changed for this and is now the default)
MC_CA02_sp <- misclass_fn(Dat.train_sp, Dat.test_sp, method="ca0", axes=2, residualised="Species")
MC_PCO2_sp <- misclass_fn(Dat.train_sp, Dat.test_sp, d=list_of_distance_matrices_sp, method="pco", axes=2, residualised="Species")
MC_CAP2_sp <- misclass_fn(Dat.train_sp, Dat.test_sp, d=list_of_distance_matrices_sp, method="cap", k=2, mp=95, axes=2, residualised="Species")

MC_CA02_CC <- misclass_fn(Dat.train_CC, Dat.test_CC, method="ca0", axes=2, residualised="CC")
MC_PCO2_CC <- misclass_fn(Dat.train_CC, Dat.test_CC, d=list_of_distance_matrices_CC, method="pco", axes=2, residualised="CC")
MC_CAP2_CC <- misclass_fn(Dat.train_CC, Dat.test_CC, d=list_of_distance_matrices_CC, method="cap", k=2, mp=95, axes=2, residualised="CC")

# amino acid level data # run 1/11/2022
MC_CA02_aa <- misclass_fn(Dat.train_aa, Dat.test_aa, method="ca0", axes=2)
MC_PCO2_aa <- misclass_fn(Dat.train_aa, Dat.test_aa, d=list_of_distance_matrices_aa, method="pco", axes=2)
MC_CAP2_aa <- misclass_fn(Dat.train_aa, Dat.test_aa, d=list_of_distance_matrices_aa, method="cap", k=2, mp=95, axes=2)

# amino acid level and CC residualisation # run 2/11/2022
MC_CA02_aa_CC <- misclass_fn(Dat.train_aa_CC, Dat.test_aa_CC, method="ca0", axes=2, residualised="CC")
MC_PCO2_aa_CC <- misclass_fn(Dat.train_aa_CC, Dat.test_aa_CC, d=list_of_distance_matrices_aa_CC, method="pco", axes=2, residualised="CC")
MC_CAP2_aa_CC <- misclass_fn(Dat.train_aa_CC, Dat.test_aa_CC, d=list_of_distance_matrices_aa_CC, method="cap", k=2, mp=95, axes=2, residualised="CC")

# Phandango distance matrix # run 19/04/23
MC_CAP2_Ph <- misclass_fn(Dat.train, Dat.test, method="cap", d=list_of_distance_matrices_Phandango, k=2, mp=95, axes=2)

## try different mp values to determine number of axes   # run 27/10/2022
mp_list <- list(0,50,85,95,99)
MC_CAP_mp <- map(mp_list, ~misclass_fn(Dat.train, Dat.test,d=list_of_distance_matrices, method="cap", k=2, mp=., axes=2))
names(MC_CAP_mp) <- map(mp_list, ~paste("CAP",.,sep="_")) |> unlist()

## Jejuni only data - check have the correct dat, Dat.train, Dat.test, and list_of_distance_matrices # run 27/10/2022
MC_CA02_j <- misclass_fn(Dat.train_j, Dat.test_j, method="ca0", axes=2)
MC_PCO2_j <- misclass_fn(Dat.train_j, Dat.test_j, d=list_of_distance_matrices, method="pco", axes=2)
MC_CAP2_j <- misclass_fn(Dat.train_j, Dat.test_j, d=list_of_distance_matrices, method="cap", k=2, mp=95, axes=2)



# merge results
MC_results <- c(list(CA01=MC_CA01, 
                   CA02=MC_CA02, 
                   PCO1=MC_PCO1, 
                   PCO2=MC_PCO2, 
                   CAP1=MC_CAP1, 
                   CAP2=MC_CAP2, 
                   CA0_CC=MC_CA02_CC,
                   PCO_CC=MC_PCO2_CC, 
                   CAP_CC=MC_CAP2_CC,
                   CA0_aa=MC_CA02_aa,
                   PCO_aa=MC_PCO2_aa,
                   CAP_aa=MC_CAP2_aa, 
                   CA0_sp=MC_CA02_sp,
                   PCO_sp=MC_PCO2_sp, 
                   CAP_sp=MC_CAP2_sp,
                   CA0_CC_aa=MC_CA02_aa_CC,
                   PCO_CC_aa=MC_PCO2_aa_CC,
                   CAP_CC_aa=MC_CAP2_aa_CC,
                   CAP_Ph=MC_CAP2_Ph,
                   CA0_j=MC_CA02_j, 
                   PCO_j=MC_PCO2_j, 
                   CAP_j=MC_CAP2_j), MC_CAP_mp)
save(MC_results, file="../CAP_Data/results/MC_results.Rdata") # run 1/11/2022, updated 2/11/22, updated 19/4/23
#load("../CAP_Data/results/MC_results.Rdata")

# now can plot - plot_results.R




## 2. Test individual tree prediction accuracy - (CAP_95 and CAP_CC_95 only) cgMLST  # run 18/11/23, updated 26/04/23
load("../CAP_Data/results/results_cgMLST.Rdata")   # from tree_data.R
# source("methods/misclassification.R") # if didn't load earlier
MC_all_trees <- misclass_tree_fn(results_cgMLST, class=Source)
#save(MC_all_trees, file="../CAP_Data/results/MC_all_trees.Rdata")  # 25/11/2022, updated 26/04/23
load(file="../CAP_Data/results/MC_all_trees.Rdata")

# now can plot ("plot_results.R")




ca <- MC_results |> map_dfr("conf.av") |> 
  mutate(method_long = names(MC_results),.before=1) |> as.data.frame() |> 
  filter(!method_long %in% c("CAP_95","CAP_CC_95"))
cs <- MC_results |> map_dfr("conf.se") |> 
  mutate(method_long = names(MC_results),.before=1) |> as.data.frame() |> 
  filter(!method_long %in% c("CAP_95","CAP_CC_95"))
a <- MC_results |> map_dfr("av") |> t() |> as.data.frame() |> rownames_to_column("method_long") |> 
  filter(!method_long %in% c("CAP_95","CAP_CC_95")) |> rename("av"="V1")
s <- MC_results |> map_dfr("se") |> t() |> as.data.frame() |> rownames_to_column("method_long")|> 
  filter(!method_long %in% c("CAP_95","CAP_CC_95"))|> rename("se"="V1")

results <- reduce(list(ca,cs,a,s), right_join, by="method_long", suffix=c(".av",".se"))|> 
  mutate(
    method = substr(method_long, 1, 3) |> factor(levels=c("CA0","PCO","CAP")),
    species = ifelse(method_long |> str_detect("j"),"j","jc"),
    axes = ifelse(method_long |> str_detect("1"),1,2),
    residualised = ifelse(method_long |> str_detect("CC"),"CC",ifelse(method_long |> str_detect("sp"), "species","no")),
    level = ifelse(method_long |> str_detect("aa"),"aa","nt"), 
    mp = ifelse(method_long |> str_extract("_[0-9]+") |> is.na(), 95, method_long |> str_extract("[0-9]+") |> as.numeric()),  
    .before=2) |> 
  pivot_longer(starts_with(c("Cattle","Chicken","Sheep")), names_to = "T_P", values_to = "p") |> rowwise() |> 
  mutate(truth = str_split(T_P, "_|\\.")[[1]][1], prediction = str_split(T_P, "_|\\.")[[1]][2])

# which had the best results?
results |> filter(truth==prediction) |> group_by(method_long, truth) |> summarise(max(p)) |> ungroup() |> rename("p"=`max(p)`) |> arrange(desc(p))
# method truth       p
# CA0    Chicken 0.873
# PCO    Chicken 0.868
# CAP    Chicken 0.868
# CAP    Sheep   0.800
# CA0    Sheep   0.800
# PCO    Sheep   0.794
# CAP    Cattle  0.626
# PCO    Cattle  0.619
# CA0    Cattle  0.619

results |> filter(truth==prediction & method=="CAP") |> group_by(method_long, truth) |> summarise(max(p)) |> ungroup() |> rename("p"=`max(p)`) |> arrange(desc(p))
# now find tree info for CAP_95 and for CAP_CC_95


# 3. Check against tree data   # run 18/11/2022
dat <- results_cgMLST |> 
  select(method, Fold, id, tree, Source, prediction) |> 
  group_by(method, id, Source) |> 
  mutate(Correct = ifelse(Source==prediction,1,0)) |> 
  select(-c(Correct,Fold)) |> mutate(prediction2 = 1) |> 
  pivot_wider(names_from = prediction, values_from = prediction2) |> 
  arrange(id) |> 
  summarise(Cattle=sum(Cattle,na.rm=TRUE), Chicken=sum(Chicken,na.rm=TRUE), Sheep=sum(Sheep,na.rm=TRUE)) |> 
  group_by(method, id, Source) |> 
  pivot_longer(cols=c(Cattle,Chicken,Sheep), names_to = "tree_prediction", values_to = "n") |> 
  summarise(prediction=tree_prediction[which.max(n)])
 
a <- dat |> mutate(True = sum(Source==prediction)) |> group_by(method, Source) |>
summarise(N=n(),n=sum(True),prop=n/N)

b <- results |> filter(method_long %in% c("CA02", "PCO2", "CAP2", "CA0_CC", "PCO_CC", "CAP_CC")) |> 
  filter(!T_P %in% c("Cattle_Cattle.se", "Chicken_Chicken.se", "Sheep_Sheep.se")) |>
  select(-c(species, axes, residualised, level, mp, T_P)) |> filter(truth==prediction) |> rename(Source = truth) |> 
  select(-method) |> rename(method=method_long) |> mutate(method=factor(method))
b$method <- recode_factor(b$method, CA02 = "CA0")
b$method <- recode_factor(b$method, PCO2 = "PCO")
b$method <- recode_factor(b$method, CAP2 = "CAP")

right_join(a, b, by=c("method", "Source")) |> select(method, Source, prop, p) |> rename(p_tree = prop, p_misscl = p)
## they should be the same but they aren't (they are almost the same!)
# method Source  p_tree p_misscl
# CA0    Cattle   0.607    0.607
# CA0    Chicken  0.844    0.844
# CA0    Sheep    0.743    0.742  *
# CA0_CC Cattle   0.585    0.584
# CA0_CC Chicken  0.862    0.863  *
# CA0_CC Sheep    0.8      0.800
# CAP    Cattle   0.607    0.619  **
# CAP    Chicken  0.844    0.844
# CAP    Sheep    0.749    0.753  **
# CAP_CC Cattle   0.612    0.612
# CAP_CC Chicken  0.867    0.868  *
# CAP_CC Sheep    0.776    0.777  *
# PCO    Cattle   0.607    0.607
# PCO    Chicken  0.849    0.853  **
# PCO    Sheep    0.733    0.732  *
# PCO_CC Cattle   0.599    0.598  *
# PCO_CC Chicken  0.867    0.868  *
# PCO_CC Sheep    0.794    0.794
