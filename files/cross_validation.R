## calculate (1) misclassification rates of the individual trees according to the presence of absent levels 
# and (2) as overal misclassification rates of the whole forests for the different methods
# Compare cgMLST and MLST

# Load libraries and functions
source("methods/libs_fns.R")
source("methods/misclassification.R")
library(caret)   # for createFolds

# Load data - jejuni and coli combined
load("../CAP_Data/data/cgMLST_dat.RData") # SACNZ cgMLST data set (jejuni and coli)
Dat_jc <- cgMLST |> mutate(across(everything(), factor)) |> 
  mutate(Source = factor(Source, levels = c("Beef", "Poultry", "Sheep", "Human"))) |>
  mutate(Source = fct_recode(as.factor(Source), Cattle="Beef", Chicken="Poultry")) |> 
  filter(Source != "Human") |> droplevels()

load("../CAP_Data/data/list_of_distance_matrices.RData") #  distance matrices of hamming distances among unique alleles for each gene
list_of_distance_matrices <- map(list_of_distance_matrices, as.matrix)

# Create data splits
set.seed(3)
flds <- createFolds(y=Dat_jc$Source, k=10) 
Dat.train <- map(flds, ~slice(Dat_jc, {-.}))  # list of k training data sets
Dat.test <- map(flds, ~slice(Dat_jc, .))    # list of k test data sets

## 1. Test individual tree prediction accuracy - cgMLST
load("../CAP_Data/results/results_cgMLST_jc.Rdata")
MC_all_tree <- misclass_tree_fn(results_cgMLST_jc)
MC_all_tree
# now can plot 


## 2. Test overall misclassification rates - cgMLST
MC_CA.zero1 <- misclass_fn(Dat.train, Dat.test, method="ca0", axes=1)
MC_CA.zero2 <- misclass_fn(Dat.train, Dat.test, method="ca0", axes=2)
MC_PCO1 <- misclass_fn(Dat.train, Dat.test, d=list_of_distance_matrices, method="pco", axes=1)
MC_PCO2 <- misclass_fn(Dat.train, Dat.test, d=list_of_distance_matrices, method="pco", axes=2)
MC_CAP1 <- misclass_fn(Dat.train, Dat.test, d=list_of_distance_matrices, method="cap", axes=1)
MC_CAP2 <- misclass_fn(Dat.train, Dat.test, d=list_of_distance_matrices, method="cap", axes=2)
MC_results <- list(CA.zero1=MC_CA.zero1, PCO1=MC_PCO1, PCO2=MC_PCO2, CAP1=MC_CAP1, CAP2=MC_CAP2)
save(MC_results, file="../CAP_Data/results/MC_results.Rdata") # run ???
load("../CAP_Data/results/MC_results.Rdata")
# now can plot 

MC_results |> map_dfc("av") |> as.data.frame()
MC_results |> map_dfc("se") |> as.data.frame()
MC_results |> map_dfr("conf.av") |> mutate(method = names(MC_results),.before=1)
MC_results |> map_dfr("conf.se") |> mutate(method = names(MC_results),.before=1)


## 3. try different mp values to determine number of axes
#mp_list <- list(0,50,85,95,99)
mp_list <- list(95)
MC_CAP <- map(mp_list, function(z) map2(Dat.train, Dat.test, ~misclass_fn(.x,.y,mp=z,d=list_of_distance_matrices, method="cap")))



# 4. Check against tree data
dat <- results_cgMLST_jc |> 
  select(method, Fold, row, tree, Source, prediction) |> 
  group_by(method, row, Source) |> 
  mutate(Correct = ifelse(Source==prediction,1,0)) |> 
  select(-c(Correct,Fold)) |> mutate(prediction2 = 1) |> 
  pivot_wider(names_from = prediction, values_from = prediction2) |> 
  arrange(row) |> 
  summarise(Beef=sum(Beef,na.rm=TRUE), Poultry=sum(Poultry,na.rm=TRUE), Sheep=sum(Sheep,na.rm=TRUE)) |> 
  group_by(method, row, Source) |> 
  pivot_longer(cols=c(Beef,Poultry,Sheep)) |> 
  summarise(max=name[which.max(value)])
dat |> summarise(True = sum(Source==max)) |> group_by(method) |> summarise(sum(True))
dat |> mutate(True = sum(Source==max)) |> group_by(method, Source) |> summarise(N=n(),n=sum(True),prop=n/N) |> arrange(Source,desc(prop))

