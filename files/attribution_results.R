#### Source attribution of the SACNZ data (jejuni and coli) using the different methods ####

# Load libraries and functions
source("methods/libs_fns.R")
source("methods/attribution.R")
source("methods/methods_plots.R")

# Load data
load("../CAP_data/data/cgMLST_dat.RData") # SACNZ cgMLST data set (jejuni and coli)
cgMLST <- cgMLST |> mutate(across(everything(), factor)) |> 
  mutate(Source = factor(Source, levels = c("Beef", "Poultry", "Sheep", "Human"))) |>
  mutate(Source = fct_recode(as.factor(Source), Cattle="Beef", Chicken="Poultry")) 

load("../CAP_data/data/list_of_distance_matrices_all.RData") #  distance matrices for animal AND human data (thus _all)
list_of_distance_matrices <- map(list_of_distance_matrices_all, as.matrix)

# Create data splits
# Training data is the source data
Dat.train <- cgMLST |> filter(Source != "Human") |> droplevels() 
# Test data is the human data to be attributed to a source
Dat.test <- cgMLST |> filter(Source == "Human") |> droplevels() 

# Run methods (no need for folds)
attribution_jc_CA.zero1 <- attribution_fn(Dat.train, Dat.test, method="ca0", axes=1)
attribution_jc_CA.zero2 <- attribution_fn(Dat.train, Dat.test, method="ca0", axes=2)
attribution_jc_PCO1 <- attribution_fn(Dat.train, Dat.test, d=list_of_distance_matrices, method="pco", axes=1)
attribution_jc_PCO2 <- attribution_fn(Dat.train, Dat.test, d=list_of_distance_matrices, method="pco", axes=2)
attribution_jc_CAP1 <- attribution_fn(Dat.train, Dat.test, d=list_of_distance_matrices, method="cap", k=2, mp=95, axes=1)
attribution_jc_CAP2 <- attribution_fn(Dat.train, Dat.test, d=list_of_distance_matrices, method="cap", k=2, mp=95, axes=2)

# Merge data sets 
attribution_jc_cgMLST <- list(
  CA.zero1 = attribution_jc_CA.zero1,
  CA.zero2 = attribution_jc_CA.zero2,
  PCO1 = attribution_jc_PCO1,
  PCO2 = attribution_jc_PCO2,
  CAP1 = attribution_jc_CAP1,
  CAP2 = attribution_jc_CAP2) |> 
  bind_rows(.id="method") |> pivot_wider(names_from = method,values_from = Freq)

attribution_jc_cgMLST

# Save data
save(attribution_jc_cgMLST, file="../CAP_Data/results/attribution_jc_cgMLST.Rdata") # run ???
load("../CAP_Data/results/attribution_jc_cgMLST.Rdata")
# now can plot 

# plot results
attribution_plot(attribution_jc_cgMLST)
