### what is the best number of PCO axes (m)?
# initially run on a random selection of 10 genes, then repeat with jejuni data and then coli data
# finally run on the full data
# and repeat using residualised distance matrices

# Load libraries and functions
source("methods/libs_fns.R")
source("methods/helpers.R")       # Functions used internally in the methods
epsilon <- sqrt(.Machine$double.eps)
source("methods/m_statistics.R")

# Function to calculate all the statistics
run <- function(dat, d, genes){
  # for each gene, run PCO and pull out propG for each m
  PCO_dat <- map(d, get_PCO_dat)
  m_stats <- data.frame(max_m=unlist(map(PCO_dat, "max_m")),nlambdas=unlist(map(PCO_dat, "nlambdas")))
  propGs <- map(PCO_dat,"propG") |> map(function(x) if(!is.null(x)) x |> mutate(m=paste0("m",m)) |> column_to_rownames("m"))
  
  # for each gene, run CA and pull out Hs
  Hs <- get_H_dat(dat, starts_with("CAMP"), "Source")
  
  # for each gene, run CAP and pull out evs (lambdaQHQs) for each m
  eigenQHQs <- map2(PCO_dat, Hs, get_eigenQHQ)
  ev1 <- map(eigenQHQs, "lambda_QHQ") |> map(function(x) map_df(x, function(y) y[1]) |> t())
  ev2 <- map(eigenQHQs, "lambda_QHQ") |> map(function(x) map_df(x, function(y) y[2]) |> t())
  
  # for each gene, run CAP leaving one allele out and repeat, then predict it, then calculate the difference and sum together
  # expected CAP scores are from the whole dataset
  CAP_scores <- map(eigenQHQs, "CAP_score")
  # returns training data, testing data, training dist, and testing dist for each gene:allele combination
  resSS.data <- subset_data(dat, genes, d)
  # calculate number of non-zero eigenvalues
  nlambdas <- map(resSS.data, ~map(.,"train.d")) |> 
    map(function(y) { 
      if(is_empty(y)) {return(NULL)}
      map_df(y, function(x){sum(eigen_decomp(dbl_center(-0.5 * x^2),symmetric = TRUE)$values > epsilon)})
    }
    ) |> 
    map(function(z) {if(is.null(z)){return(NULL)}
      min(z)})
  # restrict CAP scores to nlambda axes
  CAP_scores_reduced <- map2(CAP_scores, nlambdas, function(x,y){
    if(is.null(y)){return(NULL)}
    x[1:y]}) |> 
    map(function(x) map(x, function(y) y |> as.data.frame() |> rownames_to_column("Allele")))
  # all the data (Qo, d.gower, lambda_B, H) for each gene:allele
  LOO_dat <- map(resSS.data, ~map(.x, get_inputs))
  # estimated scores for each LOO allele
  LOO_scores <- map2(nlambdas, LOO_dat, get_score)
  # combine data from estimated scores and expected scores and take absolute values (sign is arbitrary)
  Est_Exp <- map2(LOO_scores, CAP_scores_reduced, function(x,y) {
    map2(x, y, ~right_join(.x,.y,by="Allele", suffix = c(".LOO", ".CAP")))
  }) |> map(function(x) map(x, function(y) y |>  mutate(across(starts_with("V"),abs))))
  # calculate residual sum of squares
  res.SS <- map(Est_Exp, function(x) x |> map_df(function(y) {
    Est <- y |> select(contains("LOO"))
    Exp <- y |> select(contains("CAP"))
    sum((Est - Exp)^2)}) |> as.data.frame() |> t())
  
  # merge all results together
  results <- map2_df(res.SS, pmap(list(ev1=ev1,ev2=ev2,propGs),cbind), 
                     ~ if(!is_empty(.x)) right_join(.x |> as.data.frame() |> setNames("resSS") |> rownames_to_column("m"), .y |> rownames_to_column("m"), by="m"),
                     .id = "Gene") |> mutate(across(where(is.double), round,3)) |> mutate(m=as.numeric(str_remove(m,"m")))
  
  results
}

# Load data
load("../CAP_Data/data/cgMLST_dat.RData") # SACNZ cgMLST data set (jejuni and coli)
Dat_jc <- cgMLST |> mutate(across(everything(), factor)) |> 
  mutate(Source = factor(Source, levels = c("Beef", "Poultry", "Sheep", "Human"))) |>
  mutate(Source = fct_recode(as.factor(Source), Cattle="Beef", Chicken="Poultry")) 

load("../CAP_Data/data/list_of_distance_matrices.RData") #  distance matrices of hamming distances among unique alleles for each gene
list_of_distance_matrices <- map(list_of_distance_matrices, as.matrix)

# Subset data
small_dat <- cgMLST |> select(LabID, Source, sample(2:1344,10)) |> filter(Source != "Human") |> droplevels() 
small_genes <- small_dat |> select(starts_with("CAMP")) |> colnames()
small_d <- list_of_distance_matrices[small_genes]

## and run ...
small_results <- run(small_dat, small_d, small_genes)

#plot m=1,2,3,... against
#(i) proportion of variability explained by m axes;
#(ii) delta_1^2 the first canonical correlation squared; and 
#(iii) delta_2^2
#(iv) resSS

plotdat <- small_results |> 
  pivot_longer(cols=-c(Gene, m), names_to = "Method", values_to = "Score") |> 
  mutate(Gene = factor(Gene), Method = factor(Method)) 
plotdat |> 
  filter(Method != "lambdaB") |> 
  ggplot(aes(x=m, y=Score, colour=Gene)) + geom_line() + facet_wrap(facets="Method",  scales="free_y", nrow = 1)


# repeat for jejuni only ...
details <- read.csv("../CAP_data/data/SACNZ_referencelist.csv") |> select(LabID, Species)
j_dat <- cgMLST |> left_join(details, by="LabID") |> filter(Species == "Jejuni") |> 
  select(LabID, Source, sample(2:1344,10)) |> filter(Source != "Human") |> droplevels() 
j_genes <- j_dat |> select(starts_with("CAMP")) |> colnames()
j_d <- map2(j_genes, list_of_distance_matrices[j_genes], function(x,y, data){
  alleles <- data |> pull({{x}}) |> levels()
  d <- y[alleles,alleles] 
  d}, j_dat) |> set_names(j_genes)
## and run ...
j_results <- run(j_dat, j_d, j_genes)
## and plot ...
j_results |> 
  pivot_longer(cols=-c(Gene, m), names_to = "Method", values_to = "Score") |> 
  mutate(Gene = factor(Gene), Method = factor(Method)) |> 
  filter(Method != "lambdaB") |> 
  ggplot(aes(x=m, y=Score, colour=Gene)) + geom_line() + facet_wrap(facets="Method",  scales="free_y", nrow = 1)

# and coli only ...
details <- read.csv("../CAP_data/data/SACNZ_referencelist.csv") |> select(LabID, Species)
c_dat <- cgMLST |> left_join(details, by="LabID") |> filter(Species == "Coli") |> 
  select(LabID, Source, sample(2:1344,10)) |> filter(Source != "Human") |> droplevels() 
c_genes <- c_dat |> select(starts_with("CAMP")) |> colnames()
c_d <- map2(c_genes, list_of_distance_matrices[c_genes], function(x,y, data){
  alleles <- data |> pull({{x}}) |> levels()
  d <- y[alleles,alleles] 
  d}, c_dat) |> set_names(c_genes)
## and run ...
c_results <- run(c_dat, c_d, c_genes)
## and plot ...
c_results |> 
  pivot_longer(cols=-c(Gene, m), names_to = "Method", values_to = "Score") |> 
  mutate(Gene = factor(Gene), Method = factor(Method)) |> 
  filter(Method != "lambdaB") |> 
  ggplot(aes(x=m, y=Score, colour=Gene)) + geom_line() + facet_wrap(facets="Method",  scales="free_y", nrow = 1)


# Run on all the data
dat <- cgMLST |> filter(Source != "Human") |> droplevels() 
genes <- dat |> select(starts_with("CAMP")) |> colnames()
d <- list_of_distance_matrices[genes]
## and run ...
m_statistics <- run(dat, d, genes)
#save(results, file="../CAP_Data/results/m_statistics.RData")
load("../CAP_Data/results/m_statistics.RData")


# repeat using residualised distance matrices ...
# residualised based on species
load("../CAP_Data/data/species_distance_matrices.Rdata") #  distance matrices of hamming distances among unique alleles for each gene
species_distance_matrices <- map(species_distance_matrices, as.matrix)
# load species data then subset data, remove Human cases and cases where species is unknown
details <- read.csv("../CAP_data/data/SACNZ_referencelist.csv") |> select(LabID, Species, CC)
species_dat <- cgMLST |> 
  left_join(details, by="LabID") |> select(LabID, Source, Species, CC, 2:1344) |> 
  mutate(across(everything(), factor)) |> 
  mutate(Source = factor(Source, levels = c("Beef", "Poultry", "Sheep", "Human"))) |>
  mutate(Source = fct_recode(as.factor(Source), Cattle="Beef", Chicken="Poultry")) |> 
  filter(Source != "Human") |> filter(!is.na(Species)) |> droplevels()  # change CC for Species and vice versa
species_genes <- species_dat |> select(starts_with("CAMP")) |> colnames()
species_d <- species_distance_matrices[species_genes]
## and run ...
results_species_res <- run(species_dat, species_d, species_genes)
#save(results_species_res, file="../CAP_Data/results/m_statistics_sp_res.RData")

# residualised based on clonal complex
load("../CAP_Data/data/CC_distance_matrices.Rdata") #  distance matrices of hamming distances among unique alleles for each gene
CC_distance_matrices <- map(CC_distance_matrices, as.matrix)
# load CC data then subset data, remove Human cases and cases where species is unknown
details <- read.csv("../CAP_data/data/SACNZ_referencelist.csv") |> select(LabID, Species, CC)
CC_dat <- cgMLST |> 
  left_join(details, by="LabID") |> select(LabID, Source, Species, CC, 2:1344) |> 
  mutate(across(everything(), factor)) |> 
  mutate(Source = factor(Source, levels = c("Beef", "Poultry", "Sheep", "Human"))) |>
  mutate(Source = fct_recode(as.factor(Source), Cattle="Beef", Chicken="Poultry")) |> 
  filter(Source != "Human") |> filter(!is.na(CC)) |> droplevels()
CC_genes <- CC_dat |> select(starts_with("CAMP")) |> colnames()
CC_d <- CC_distance_matrices[CC_genes]
## and run ...
results_CC_res <- run(CC_dat, CC_d, CC_genes)
#save(results_CC_res, file="../CAP_Data/results/m_statistics_CC_res.RData")


# did residualising the distance matrices help?
m_res <- right_join(results, results_species_res, by=c("Gene", "m"), suffix = c("", "_sp.res"))
m_res <- right_join(results, results_CC_res, by=c("Gene", "m"), suffix = c("", "_CC.res"))
m_res |> filter(m %in% c(1,2)) |> mutate(across(m, factor)) |> select(Gene, m, starts_with("Var")) |> 
  pivot_longer(cols = starts_with("Var"), names_to = "Method", values_to = "VarExp") |> 
  ggplot(aes(x=m, y=VarExp, fill=Method)) + geom_boxplot() + theme_bw() + scale_fill_brewer(palette="Dark2")
m_res |> select(Gene, m, starts_with("ev"))
m_res |> filter(m %in% c(1,2)) |> select(Gene, m, starts_with("res"))


