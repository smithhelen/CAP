## Extending to two genes

# load libraries and functions
source("methods/libs_fns.R")
source("methods/tree_predictions.R")  # to pull out individual tree predictions
library(boot) # for inv.logit function

# set seed
set.seed(18)

## generate data
gen_data_2genes <- function(beta){

  # gene 1 has 10 alleles:  
  # generate PCO1 and PCO2 scores - the direction of greatest variation is along PCO1
  Q1 <- data.frame(PCO1=runif(10, -5, 5), PCO2=runif(10, -1, 1)) |> arrange(PCO2)  # arrange so that the lowest numbers prefer source 1
  # generate probability scores - the probability that an allele will belong to source one.
  p1 <- inv.logit(beta * Q1$PCO2) 
  # generate counts for each allele from 14 individuals (size=14)
  X1 <- data.frame(n1=(map(p1, ~rbinom(n=1, size=14, p=.)) |> unlist())) |> mutate(n2=14-n1)
  
  # gene 2 has 7 alleles:
  # generate PCO1 and PCO2 scores - the direction of greatest variation is along PCO1
  Q2 <- data.frame(PCO1=runif(7, -5, 5), PCO2=runif(7, -1, 1)) |> arrange(PCO2) # arrange so that the lowest numbers prefer source 1
  # generate probability scores - the probability that an allele will belong to source one.
  p2 <- inv.logit(beta * Q2$PCO2) 
  # generate counts for each allele from 20 individuals (size=20)
  X2 <- data.frame(n1=(map(p2, ~rbinom(n=1, size=20, p=.)) |> unlist())) |> mutate(n2=20-n1)
  
  # calculate distances
  d1 <- Q1 |> dist(diag = TRUE) |> as.matrix()
  d2 <- Q2 |> dist(diag = TRUE) |> as.matrix()
  
  # fill out individual row data ie what the original data would actually look like
  df <- data.frame(gene = rep(factor(1:2), c(10,7)), 
                   allele = c(factor(1:10), factor(1:7)), 
                   PCO1 = c(Q1$PCO1, Q2$PCO1), 
                   PCO2 = c(Q1$PCO2, Q2$PCO2), 
                   s1 = c(X1$n1,X2$n1), 
                   s2 = c(X1$n2,X2$n2),
                   observed = "N")
  ob <- c(sample(10,7),sample(7,4)+10)
  df[ob,"observed"] <- "Y"
  
  dat <- df |> pivot_longer(cols = starts_with("s"), names_to = "source", values_to = "n", names_prefix = "count_") |> 
    uncount(n)|> mutate(across(!where(is.double), as.factor))
  dat[dat$gene==2,"source"] <- dat[dat$gene==1,"source"]
  # check all combinations are in data
  #table(dat[dat$gene==1,]$observed, dat[dat$gene==2,]$observed
  
  dat <- dat |> filter(gene == 1) |> select(allele, observed) |> rename(gene1 = allele, observed1 = observed) |> 
    bind_cols(dat |> filter(gene == 2) |> select(allele, observed, source) |> rename(gene2 = allele, observed2 = observed)) |> 
    mutate(observed = ifelse(observed1=="Y"&observed2=="Y","Y","N")) |> select(gene1, gene2, source, observed)|> 
    mutate(across(!where(is.double), as.factor)) |> rownames_to_column("id")
  
  simdat <- list(df=dat, d=list(d1=d1,d2=d2))
  simdat
}

# the training data is a sample all of the individuals which have 'observed' alleles
# the testing data has a mix of observed and unobserved alleles
split_dat <- function(dat){
  split.id <- dat |> filter(observed == "Y") |> slice_sample(prop=0.8) |> pull(id)
  split.id
}

## encode factors
encode <- function(dat.train, dat.test, d=NULL, axes=2, mp=NULL, m=2, k=2, id=id, class, var_id){
  train.ca0 <- prepare_training_ca0(dat.train, starts_with(var_id), class=class, axes=axes)
  test.ca0 <- prepare_test_ca0(dat.test, train.ca0$extra, id=id)
  train.pco <- prepare_training_pco(dat.train, starts_with(var_id), class=class, d, axes=axes)
  test.pco <- prepare_test_pco(dat.test, train.pco$extra, id=id)
  train.cap <- prepare_training_cap(dat.train, starts_with(var_id), class=class, d, axes=axes, k=k, m=m, mp=mp)
  test.cap <- prepare_test_cap(dat.test, train.cap$extra, id=id)
  out <- list(ca0 = list(train=train.ca0, test=test.ca0),
              pco = list(train=train.pco, test=test.pco),
              cap = list(train=train.cap, test=test.cap))
  out
}

## run ranger
run_ranger <- function(dat, ntrees=500, d=d, k=2, m=2, axes=1, class="source", id="id", var_id="gene"){
  # prepare data - add id and split into training and test sets
  # the training data is a sample all of the individuals which have 'observed' alleles
  # the testing data has a mix of observed and unobserved alleles
  split.id <- split_dat(dat)
  dat.train <- dat |> filter(id %in% split.id)
  dat.test <- dat |> anti_join(dat.train)
  
  # set number of trees
  ntrees = ntrees
  
  # encode factors
  dat.scores <- encode(dat.train, dat.test, d=d, k=2, m=2, axes=1, class=class, id=id, var_id=var_id, ntrees=ntrees)
  genes <- colnames(dat |> select(starts_with("gene")))
  names(dat.scores) <- genes
  
  # find uniques
  extras <- dat.scores |> map(c("train","extra"))
  uniques <- extras |> pluck(1) |> is_unique(data=dat.test)
  
  # train model
  training <- dat.scores |> map(c("train", "training")) |> map(function(x) x |> select(-all_of(class)))
  classes <- dat.train |> pull(class)
  rf_mods <- map(training, ~ranger(classes ~ ., data=.x, oob.error = TRUE, num.trees=ntrees, respect.unordered.factors = TRUE))
  rf_mods_pco <- ranger(classes ~ ., data=dat.train |> select(PCO1), oob.error = FALSE, num.trees=ntrees, respect.unordered.factors = TRUE) # ranger with scores (not dist) for PCO method
  
  # individual tree predictions
  testing <- dat.scores |> map("test")
  tree_preds <- map2(rf_mods, testing, ~predict_by_tree(.x, .y, uniques, id=id))
  tree_preds_pco <- predict_by_tree(rf_mods_pco, dat.test |> select(all_of(id), PCO1), uniques |> rename(PCO1=allele.V1), id=id)
  answer <- tree_preds |> map(function(x) x |> left_join(dat.test |> select(all_of(id), all_of(class))))
  answer_pco <- tree_preds_pco  |> left_join(dat.test |> select(all_of(id), all_of(class)))
  #(sim_PCO$prediction == sim_PCOb$prediction) |> table() #not exactly the same
  
  list(tree_preds = answer, tree_preds_pco = answer_pco, dat.scores=dat.scores)
}

## prepare data
prepare_dat <- function(dat, d=d, k=2, m=2, axes=1, class="source", id="id"){

  # identify genes
  
  # encode factors of each gene (will be a list)
  dat.scores <- encode(dat.train, dat.test, d=list(y), k=2, m=2, axes=1, class=class, id=id, var_id=starts_with("gene")))
  names(dat.scores) <- genes
  dat.scores
}

## run ranger
run_ranger <- function(dat.scores, ntrees=500, class="source", id="id"){
  # set number of trees
  ntrees = ntrees
  
  # find uniques
  extras <- dat.scores |> map(c("train","extra"))
  uniques <- extras |> pluck(1) |> is_unique(data=dat.test)
  
  # train model
  training <- dat.scores |> map(c("train", "training")) |> map(function(x) x |> select(-all_of(class)))
  classes <- dat.train |> pull(class)
  rf_mods <- map(training, ~ranger(classes ~ ., data=.x, oob.error = TRUE, num.trees=ntrees, respect.unordered.factors = TRUE))
  rf_mods_pco <- ranger(classes ~ ., data=dat.train |> select(PCO1), oob.error = FALSE, num.trees=ntrees, respect.unordered.factors = TRUE) # ranger with scores (not dist) for PCO method
  
  # individual tree predictions
  testing <- dat.scores |> map("test")
  tree_preds <- map2(rf_mods, testing, ~predict_by_tree(.x, .y, uniques, id=id))
  tree_preds_pco <- predict_by_tree(rf_mods_pco, dat.test |> select(all_of(id), PCO1), uniques |> rename(PCO1=allele.V1), id=id)
  answer <- tree_preds |> map(function(x) x |> left_join(dat.test |> select(all_of(id), all_of(class))))
  answer_pco <- tree_preds_pco  |> left_join(dat.test |> select(all_of(id), all_of(class)))
  #(sim_PCO$prediction == sim_PCOb$prediction) |> table() #not exactly the same
  
  list(tree_preds = answer, tree_preds_pco = answer_pco, dat.scores=dat.scores)
}

mc_trees <- function(dat, class){
  df <- dat |> 
    mutate(uses_unique = if_else(uses_unique == 0, "No", "Yes")) |>
    mutate(correct = across(all_of(class))==prediction) |> 
    group_by(method,uses_unique, correct) |>
    summarise(n=n()) |>
    group_by(method,uses_unique) |>
    mutate(N = sum(n)) |>
    mutate(prop = n/N) |>     
    filter(correct==TRUE) |> 
    mutate(mc = 1-prop)
  df
}

# put it into a function
sim_fn <- function(beta, ntrees){
  # generate data
  simdat <- gen_data_2genes(beta)
  
  # run ranger predictions
  preds <- run_ranger(dat = simdat$df, ntrees, d=simdat$d, k=2, m=2, axes=1, class="source", id="id", var_id="allele")
  
  #misclassification
  mcdat <- c(preds$tree_preds, pco2=list(preds$tree_preds_pco)) |> bind_rows(.id = "method")
  mc <- mc_trees(mcdat, class="source")
  
  # store misclassification rates for observed alleles and unobserved alleles
  # output results
  out <- list(mc=mc, mc.extras=list(dat=dat, tree_preds=c(preds$tree_preds, pco2=list(preds$tree_preds_pco)), dat.scores=preds$dat.scores))
  out
}


