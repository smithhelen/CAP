### calculate variable importance values of the first 4 methods
# note this is not the best way to determine the importance of variables as it uses the OOB data 
# but is a simple way to check if the same data is being used for the different methods (not CAP yet)

# don't worry about the folds as not calculating misclassification rates.

# Load libraries and functions
source("methods/libs_fns.R")

# Load data
load("../CAP_Data/data/cgMLST_dat.RData") # SACNZ cgMLST data set (jejuni and coli)
Dat_jc <- cgMLST |> mutate(across(everything(), factor)) |> 
  mutate(Source = factor(Source, levels = c("Beef", "Poultry", "Sheep", "Human"))) |>
  mutate(Source = fct_recode(as.factor(Source), Cattle="Beef", Chicken="Poultry")) 

load("../CAP_Data/data/list_of_distance_matrices_all.RData") #  distance matrices of hamming distances among unique alleles for each gene
list_of_distance_matrices <- map(list_of_distance_matrices_all, as.matrix)

# Create data splits (no need for test data)
Dat.train <- Dat_jc |> filter(Source != "Human") |> droplevels() 

# Run methods (no need for folds)
importance_CA.zero1 <- var_impt_fn(Dat.train, method="ca0", axes=1)
importance_CA.zero2 <- var_impt_fn(Dat.train, method="ca0", axes=2)
importance_PCO1 <- var_impt_fn(Dat.train, d=list_of_distance_matrices, method="pco", axes=1)
importance_PCO2 <- var_impt_fn(Dat.train, d=list_of_distance_matrices, method="pco", axes=2) 
importance_CAP1 <- var_impt_fn(Dat.train, d=list_of_distance_matrices, method="cap", k=2, m=2, axes=1)
importance_CAP2 <- var_impt_fn(Dat.train, d=list_of_distance_matrices, method="cap", k=2, m=2, axes=2) 

## Merge data sets
importance <- bind_rows(CA.zero1 = importance_CA.zero1 %>% slice(1:10) %>% select(gene) %>% mutate(imp = c(1:10)),
                        CA.zero2 = importance_CA.zero2 %>% slice(1:10) %>% select(gene) %>% mutate(imp = c(1:10)),
                        PCO1 = importance_PCO1 %>% slice(1:10) %>% mutate(gene = gene %>% substr(0,8)) %>% select(gene) %>% mutate(imp = c(1:10)),
                        PCO2 = importance_PCO2 %>% slice(1:10) %>% mutate(gene = gene %>% substr(0,8)) %>% select(gene) %>% mutate(imp = c(1:10)),
                        CAP1 = importance_CAP1 %>% slice(1:10) %>% mutate(gene = gene %>% substr(0,8)) %>% select(gene) %>% mutate(imp = c(1:10)),
                        CAP2 = importance_CAP2 %>% slice(1:10) %>% mutate(gene = gene %>% substr(0,8)) %>% select(gene) %>% mutate(imp = c(1:10)),
                        .id="method") %>% pivot_wider(names_from = method, values_from = imp) %>% as.data.frame()
importance
