# Create data files to test against PRIMERs CAP routine
# libraries
library(tidyverse)
library(ranger)
library(FactoMineR)
library(factoextra)

# load data
load("../CAP_data/PRIMER/Dat_j.RData") # SACNZ cgMLST data set for jejuni only
load("../CAP_data/PRIMER/list_of_distance_matrices_j.RData") #  distance matrices of hamming distances among unique alleles for each gene
Dat_j <- Dat_j |> mutate(across(everything(), factor)) 
list_of_distance_matrices_j <- map(list_of_distance_matrices_j, as.matrix)

# create small set of test data
vars <- colnames(Dat_j)[c(2,5,14,20)]
data <- Dat_j |> select(LabID, Source, all_of(vars))
dist <- list_of_distance_matrices_j[vars]

# create data splits
set.seed(2)
split <- round(runif(round(0.2*nrow(data),0),1,nrow(data)),0)
Dat.train <- slice(data, -split)  
Dat.test <- slice(data, split)    

# calculate capX scores for CAP routine in PRIMER
capX.fun <- function(genes, data, source) {
  gene <- data |> pull({{genes}}) |> droplevels()
  Source <- data |> pull({{source}})
  ct <- table(Allele=gene,Source)
  ca1 <- CA(ct, ncp=2, graph=FALSE)
  capX <- get_ca_row(ca1)$coord |> as.data.frame() |> add_column(NA) |> mutate(Allele = levels(gene), Group = "Train")
  colnames(capX)[3] <- ""
  write.csv(capX, file=paste0("PRIMER_verify/capX_train_",genes,".csv"), na="")
}
#map(vars, capX.fun, data=Dat.train, Source)

capd.fun <- function(d, gene, train, test){
  group <- rep("Train", times=ncol(d))
  order <- sort(colnames(d))
  group[order %in% setdiff(test |> pull({{gene}}), train |> pull({{gene}}))] <- "Test"
  dist <- d[order,order] |> as.data.frame() |> rbind("") |> rbind(Allele = order, Group = group)
  rownames(dist)[nrow(dist)-2] <- ""
  write.csv(dist, file=paste0("PRIMER_verify/d_",{{gene}},".csv"), na="")
}
#map2(dist, names(dist), ~capd.fun(.x,.y, Dat.train, Dat.test))

# import into PRIMER and run program. 
# See "../CAP_data/PRIMER/PRIMER_guide_to_running_CAP.pdf" for details on running the CAP and
# "../CAP_data/PRIMER/PRIMER_output_results.pdf" for how to output the results.

# adjust R CAP method to give output we need to check
source("methods/cap.R") 

P_prepare_test_cap <- function(data, extra, id) {
  id <- data |> select({{id}})
  var_cols <- data |> select(any_of(names(extra)))
  newdata_score <- map2(var_cols, extra, P_impute_score_cap)
  # extra line for comparing Primer results
  CAP_score <- map(newdata_score,"var_level_score")
  output <- map(newdata_score,"test_score")
  newdata_pred <- map2_dfc(output, names(output), ~ .x |> set_names(paste(.y, names(.x), sep=".")))
  newdata_pred <- bind_cols(id, newdata_pred)
  # list to output extra results for Primer compare
  list(newdata_pred=newdata_pred, CAP_score=CAP_score)
}

P_impute_score_cap <- function(var, extra) {
  var_levels <- pluck(extra, "var_levels")
  var <- droplevels(var)
  new.var_levels <- setdiff(levels(var), var_levels)
  d <- pluck(extra, "d")
  diag.B.train <- pluck(extra, "diag.B.train")
  lambda_B <- pluck(extra, "lambda_B")
  Qo <- pluck(extra, "Qo")
  C_score <- pluck(extra, "C_score")
  U <- pluck(extra, "U")
  lambda_QHQ <- pluck(extra, "lambda_QHQ")
  new_scores <- map_df(new.var_levels, ~predict_cap(new.var_level = {.}, d, diag.B.train, Qo, lambda_B, lambda_QHQ, U))
  # extra line for comparing Primer results
  var_level_score <- rbind(data.frame(C_score) |> rownames_to_column("Var_level"), new_scores)
  # extra added to list to output extra results for Primer compare
  list(var_level_score = var_level_score,
       test_score = data.frame(Var_level = var) |> 
         left_join(var_level_score, by = "Var_level") |>
         select(-Var_level))
}

prepped_train <- prepare_training_cap(Dat.train, starts_with("CAMP"), "Source", dist, k=2, m=2, axes=2)
prepped_test <- P_prepare_test_cap(Dat.test, prepped_train$extra, "LabID")
R_results <- prepped_test$CAP_score

# compare results against CAP_CAMPxxxx.txt files
read_PRIMER_txt <- function(gene){
  data <- read_tsv(paste0("../CAP_data/PRIMER/CAP_",{{gene}},".txt")) |>
    select(c(starts_with("CAP")),Allele) |> as.data.frame() |>
    rename_with(~gsub("CAP","V",.x,fixed=TRUE), starts_with("CAP")) |>
    relocate(Allele)
  data
}

P_results <- map(vars, read_PRIMER_txt) |> setNames(vars)

check_diff <- function(R_result, P_result){
  vars <- as.list(colnames(R_result |> select(starts_with("V"))))
  check <- function(var){
    R <- R_result |> select({{var}})
    P <- P_result |> select({{var}})
    sum(pmin(abs(R + P),abs(R-P))) |> set_names({{var}})
  } 
  diff <- map(vars, check)
  diff |> unlist()
}
map2(P_results, R_results, check_diff)


