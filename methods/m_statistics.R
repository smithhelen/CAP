### what is the best number of PCO axes (m)?
## all the functions

# for each gene, run PCO and pull out propG for each m
get_PCO_dat <- function(d.train){
  A.train <- -0.5 * d.train^2
  B.train <- dbl_center(A.train)  #B is same as Gowers matrix G
  eigen_B <- eigen_decomp(B.train, symmetric=TRUE)
  Qo <- eigen_B$vectors
  lambdaB <- eigen_B$values
  m <- 1:length(lambdaB)
  #number of positive/nonzero eigenvalues
  nlambdas <- sum(lambdaB > epsilon)
  if (nlambdas == 0) {
    # No non-zero eigenvectors
    return(NULL)
  }
  #propG proportion of variation explained by m axes
  VarExp <- cumsum(lambdaB/sum(lambdaB)*100) |> round(6)
  #number of eigenvalues which don't go over 100%
  if(nlambdas != max(m)){
    max_m <- max(which(VarExp<=100))
  }
  else max_m <- m
  propG <- lambdaB |> as.data.frame() |> rownames_to_column("m") |> cbind(VarExp)
  list(Qo=Qo,
       propG = propG, 
       nlambdas = nlambdas, 
       max_m=max_m)
}

# for each gene, run CA and pull out Hs
get_H <- function(var, class){
  var_levels <- droplevels(var)
  if(var_levels |> unique() |> length() < 2){
    return(NULL)
  }
  ct <- table(Var_level=var_levels,Class=class)
  H <- hat(ct, k=2)
  H
}
get_H_dat <- function(data, vars, class){
  var_cols <- data |> select({{vars}})
  classes   <- data |> pull({{class}})
  H <- map(var_cols, get_H, class = classes)
  H
}

# for each gene, run CAP and pull out evs (lambdaQHQs) for each m
get_eigenQHQ <- function(Q, H){
  if(is.null(Q)){return(NULL)}
  Qo <- Q |> pluck("Qo")
  m <- nrow(Qo)
  #m <- Q |> pluck("max_m")
  Qm <- map(1:m, function(x) Qo[,1:x] ) |> set_names(paste0("m",1:m) )
  QHQ <- map(Qm, function(x) t(x) %*% H %*% x)
  eigenQHQs <- map(QHQ, function(x) {eigen_decomp(x, symmetric=TRUE)})
  lambda_QHQ <- map(eigenQHQs, "values")
  U <- map(eigenQHQs, "vectors")
  CAP_score <- map2(Qm, U, function(x,y) {
    c <- x %*% y
    rownames(c) <- rownames(Qo)
    return(c)
  })
  list(lambda_QHQ=lambda_QHQ, 
       CAP_score=CAP_score)
}

# resSS
# for each gene, run CAP leaving one allele out and repeat, then predict it, then calculate the difference and sum together
# create training set which is all alleles bar one and test set which is the LOO allele
filter_alleles <- function(a, dat, d) {
  train.dat <- dat |> filter(if_all(starts_with("CAMP"), ~. != {{a}})) |> droplevels()
  test.dat <- dat  |> filter(if_all(starts_with("CAMP"), ~. == {{a}})) |> droplevels()
  train.d <- d[colnames(d) != {{a}}, colnames(d) != {{a}}]
  test.d <- d[{{a}}, colnames(d) != {{a}}]
  list(train=train.dat, test=test.dat, train.d=train.d, test.d=test.d)
}
# pull out named list of alleles
subset_alleles <- function(dat, dist) {
  if(is.null(dim(dist))){return(NULL)}
  d <- dist
  alleles <- dat |> select(starts_with("CAMP")) |> pull() |> levels() |> as.list() 
  names(alleles) <- alleles
  map(alleles, filter_alleles, dat, d)
}
# pull out individual genes as dataframes
subset_data <- function(dat, genes, dist) {
  single.data <- map(genes, ~ dat |> select(LabID, Source, all_of(.x))) |> setNames(genes)
  data.sets <- map2(single.data, dist, subset_alleles)
}
# get input arguments
get_inputs <- function(allele) {
  train <- allele$train
  #cat(colnames(train)[3])
  d.train <- allele$train.d
  test <- allele$test
  d.test <- allele$test.d
  #pco bit
  A.train <- -0.5 * d.train^2
  B.train <- dbl_center(A.train)  #B is same as Gowers matrix G
  eigen_B <- eigen_decomp(B.train, symmetric=TRUE)
  lambda_B <- eigen_B$values
  nlambdas <- sum(lambda_B > epsilon)
  if (nlambdas == 0) {
    return(NULL)
  }
  Qo <- eigen_B$vectors
  d.gower <- diag(B.train) - (d.test ^ 2)
  #ca bit
  ct <- train |>select(3,2) |> table()
  H <- hat(ct, k=2)
  #repeat for each m
  list(Qo=Qo,
       d.gower = d.gower, 
       lambda_B = lambda_B, 
       H=H)
}
# estimate CAP scores for the LOO allele
est_cap <- function(m, Qo, d.gower, lambda_B, H){
  lambda_Bm <- map(lambda_B, ~.x[1:m])
  Qm <- map(Qo, ~.x[,1:m])
  QHQ <- map2(Qm, H, ~t(.x) %*% .y %*% .x)
  eigenQHQs <- map(QHQ, ~eigen_decomp(.x, symmetric=TRUE))
  U <- map(eigenQHQs, "vectors")
  newQo <- pmap(list(d.gower, Qm, lambda_Bm), ~ ({..1} %*% {..2}) / (2*{..3}))
  newCAPscore <- map2(newQo, U, ~ .x %*% .y) 
  map_df(newCAPscore, as.data.frame, .id="Allele")
}
# get CAP scores for the training data
get_score <- function(nlambda, data){
  if(is.null(nlambda)){return(NULL)}
  Qo <- data |> map("Qo")
  d.gower <- data |> map("d.gower")
  lambda_B <- data |> map("lambda_B")
  H <- data |> map("H")
  map(seq_len(nlambda), ~est_cap(m=.x, Qo, d.gower, lambda_B, H) ) |> set_names(paste0("m",seq_len(nlambda)))
}

