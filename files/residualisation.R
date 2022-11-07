## Load libraries and functions ##
source("methods/libs_fns.R")
source("methods/residualise.R")
library(fastDummies)

## Load data ##
# allelic information
load("../CAP_Data/data/cgMLST_dat.RData") # SACNZ cgMLST data set (jejuni and coli)

# distance information (Hamming)
load("../CAP_Data/data/list_of_distance_matrices.RData") #  distance matrices of hamming distances among unique alleles for each gene
list_of_distance_matrices <- map(list_of_distance_matrices, as.matrix)


## 1. for residualisation based on species ## ##################################################################################
# load species data
species <- read.csv("../CAP_data/data/SACNZ_referencelist.csv") |> select(LabID, Species)

# subset data, remove Human cases and cases where species is unknown
Dat_species <- cgMLST |> 
  left_join(species, by="LabID") |> select(LabID, Source, Species, 2:1344) |> 
  mutate(across(everything(), factor)) |> 
  mutate(Source = factor(Source, levels = c("Beef", "Poultry", "Sheep", "Human"))) |>
  mutate(Source = fct_recode(as.factor(Source), Cattle="Beef", Chicken="Poultry")) |> 
  filter(Source != "Human") |> filter(!is.na(Species)) |> droplevels()
Dat_species |> as_tibble()

# model, hat, and identity matrices
X <- model.matrix(~Dat_species$Species)    # because there are only 2 groups we don't need to worry about adjusting for sample sizes
H <- X %*% solve(t(X) %*% X) %*% t(X)
I <- H |> ncol() |> diag()

# residualise each distance matrix
var_cols <- Dat_species |> select(starts_with("CAMP"))
species_distance_matrices <- map2(var_cols, list_of_distance_matrices, residualise, I, H)
#save(species_distance_matrices, file="../CAP_Data/data/species_distance_matrices.Rdata")
load("../CAP_Data/data/species_distance_matrices.Rdata") #  distance matrices of hamming distances among unique alleles for each gene
species_distance_matrices <- map(species_distance_matrices, as.matrix)

# run CAP cross-validation (specifying species in the model) (cross_validation.R)


## 2. for residualisation based on CC ## #####################################################################################
# load clonal complex data
CC <- read.csv("../CAP_data/data/SACNZ_referencelist.csv") |> select(LabID, CC)

# subset data, remove Human cases and cases where CC is unknown
Dat_CC <- cgMLST |> 
  left_join(CC, by="LabID") |> select(LabID, Source, CC, 2:1344) |> 
  mutate(across(everything(), factor)) |> 
  mutate(Source = factor(Source, levels = c("Beef", "Poultry", "Sheep", "Human"))) |>
  mutate(Source = fct_recode(as.factor(Source), Cattle="Beef", Chicken="Poultry")) |> 
  filter(Source != "Human") |> filter(!is.na(CC)) |> droplevels() |> arrange("LabID")
Dat_CC |> as_tibble()

## model, hat, and identity matrices
# calculate relative frequencies
p <- Dat_CC |> group_by(CC) |> summarise(p=1/n())

# convert CC variable to dummy
CC_dummy <- dummy_cols(Dat_CC |> select(LabID,Source,CC), select_columns="CC") |> arrange(CC) 

# substitute relative frequency values with dummy 1s
contrast_fn <- function(dat, i){
  data.frame(dat[,i]-dat[,i+1]) |> set_names(colnames(dat)[i])}
CC_dat <- CC_dummy |> select(-c(LabID, Source, CC)) * rep(p$p, each=nrow(Dat_CC))
CC_contrasts <- cbind(
  LabID=CC_dummy$LabID,
  map_dfc(1:(nrow(p)-1), ~contrast_fn(CC_dat, i=.))
)

# align rows of X with rows of B
X <- CC_contrasts |> arrange(LabID) |> select(-LabID) |> as.matrix() 
H <- X %*% solve(t(X) %*% X) %*% t(X)
I <- H |> ncol() |> diag()

# residualise each distance matrix
var_cols <- Dat_CC |> select(starts_with("CAMP"))
CC_distance_matrices <- map2(var_cols, list_of_distance_matrices, residualise, I, H)
#save(CC_distance_matrices, file="../CAP_Data/data/CC_distance_matrices.Rdata")
#load("../CAP_Data/data/CC_distance_matrices.Rdata") #  distance matrices of hamming distances among unique alleles for each gene
CC_distance_matrices <- map(CC_distance_matrices, as.matrix)


## 3. amino acid data with residualisation based on CC ## #####################################################################################
# load amino acid alleleic and distance information
load("../CAP_data/data/aa_allele_dat.RData")
load("../CAP_data/data/aa_distance_matrices.RData")

# subset data, remove Human cases and cases where CC is unknown
Dat_aa_CC <- aa_allele_dat |> 
  left_join(CC, by="LabID") |> select(LabID, Source, CC, starts_with("CAMP")) |> 
  mutate(across(everything(), factor)) |> 
  mutate(Source = factor(Source, levels = c("Beef", "Poultry", "Sheep", "Human"))) |>
  mutate(Source = fct_recode(as.factor(Source), Cattle="Beef", Chicken="Poultry")) |> 
  filter(Source != "Human") |> filter(!is.na(CC)) |> droplevels() |> arrange("LabID")
Dat_aa_CC |> as_tibble()

## model, hat, and identity matrices
# calculate relative frequencies
p <- Dat_aa_CC |> group_by(CC) |> summarise(p=1/n())

# convert CC variable to dummy
CC_dummy <- dummy_cols(Dat_aa_CC |> select(LabID,Source,CC), select_columns="CC") |> arrange(CC) 

# substitute relative frequency values with dummy 1s
contrast_fn <- function(dat, i){
  data.frame(dat[,i]-dat[,i+1]) |> set_names(colnames(dat)[i])}
CC_dat <- CC_dummy |> select(-c(LabID, Source, CC)) * rep(p$p, each=nrow(Dat_aa_CC))
CC_contrasts <- cbind(
  LabID=CC_dummy$LabID,
  map_dfc(1:(nrow(p)-1), ~contrast_fn(CC_dat, i=.))
)

# align rows of X with rows of B
X <- CC_contrasts |> arrange(LabID) |> select(-LabID) |> as.matrix() 
H <- X %*% solve(t(X) %*% X) %*% t(X)
I <- H |> ncol() |> diag()

# residualise each distance matrix
var_cols <- Dat_aa_CC |> select(starts_with("CAMP"))
CC_distance_matrices_aa <- map2(var_cols, aa_distance_matrices, residualise, I, H)
CC_distance_matrices_aa <- map(CC_distance_matrices_aa, as.matrix)
#save(CC_distance_matrices_aa, file="../CAP_Data/data/CC_distance_matrices_aa.Rdata")  #run 2/11/2022
#load("../CAP_Data/data/CC_distance_matrices_aa.Rdata") #  distance matrices of hamming distances among unique alleles for each gene


#### run CAP cross-validation (specifying CC in the model) (cross_validation.R) ####

