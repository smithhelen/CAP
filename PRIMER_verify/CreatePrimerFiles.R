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
prepare_test_cap <- function(data, extra) {
  gene_cols <- data |> select(any_of(names(extra)))
  newdata_score <- map2(gene_cols, extra, impute_score_cap)
  CAP_score <- map(newdata_score,"allele_score")
  output <- map(newdata_score,"test_score")
  newdata_pred <- map2_dfc(output, names(output), ~ .x |> set_names(paste(.y, names(.x), sep=".")))
  list(newdata_pred=newdata_pred, CAP_score=CAP_score)
}
impute_score_cap <- function(gene, extra) {
  alleles <- pluck(extra, "alleles")
  gene <- droplevels(gene)
  new.alleles <- setdiff(levels(gene), alleles)
  d <- pluck(extra, "d")
  B.train <- pluck(extra, "B.train")
  lambda_B <- pluck(extra, "lambda_B")
  Q <- pluck(extra, "Q")
  Q_score <- pluck(extra, "Q_score")
  U <- pluck(extra, "U")
  lambda_QHQ <- pluck(extra, "lambda_QHQ")
  axes <- pluck(extra, "axes")
  new_scores <- map_df(new.alleles, ~predict_cap(new.allele = {.}, d, B.train, Q, lambda_B, lambda_QHQ, U, axes))
  allele_score <- rbind(data.frame(Q_score) |> rownames_to_column("Allele"), new_scores)
  list(allele_score = allele_score,
       test_score = data.frame(Allele = gene) |> 
         left_join(allele_score, by = "Allele") |>
         select(-Allele))
}
predict_cap <- function(new.allele, d, B.train, Q, lambda_B, lambda_QHQ, U, axes) {
  # new.allele is length one
  d.new <- d[new.allele, colnames(B.train)]
  d.gower <- diag(B.train) - (d.new ^ 2)
  X <- sweep(Q, 2, sqrt(abs(lambda_B[1:axes])), "*")
  newX <- d.gower %*% X / (2 * lambda_B[1:axes])
  newQ <- newX / sqrt(abs(lambda_B[1:axes]))
  newQscore <- newQ %*% U * sqrt(abs(lambda_QHQ))
  new_allele_score <- data.frame(Allele = new.allele, newQscore)
  new_allele_score
}

prepped_train <- prepare_training_cap(Dat.train, starts_with("CAMP"), Source, dist, axes=2)
prepped_test <- prepare_test_cap(Dat.test, prepped_train$extra)
R_results <- prepped_test$CAP_score

# compare results against CAP_CAMPxxxx.txt files
read_PRIMER_txt <- function(gene){
  data <- read_tsv(paste0("PRIMER_verify/CAP_",{{gene}},".txt")) |>
    select(!(starts_with("X") | Group)) |> as.data.frame() |>
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


