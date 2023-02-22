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
  
  simdat <- list(df=dat, d=list(gene1=d1,gene2=d2))
  simdat
}

# the training data is a sample all of the individuals which have 'observed' alleles
# the testing data has a mix of observed and unobserved alleles
split_dat <- function(dat){
  split.id <- dat |> filter(observed == "Y") |> slice_sample(prop=0.8) |> pull(id)
  split.id
}

## encode factors
encode <- function(Dat.train, Dat.test, d=NULL, axes=2, mp=NULL, m=2, k=2, ntrees=500, id=id, class, var_id){
  train.ca0 <- prepare_training_ca0(Dat.train, starts_with(var_id), class=class, axes=axes)
  test.ca0 <- prepare_test_ca0(Dat.test, train.ca0$extra, id=id)
  train.pco <- prepare_training_pco(Dat.train, starts_with(var_id), class=class, d, axes=axes)
  test.pco <- prepare_test_pco(Dat.test, train.pco$extra, id=id)
  train.cap <- prepare_training_cap(Dat.train, starts_with(var_id), class=class, d, axes=axes, k=k, m=m, mp=mp)
  test.cap <- prepare_test_cap(Dat.test, train.cap$extra, id=id)
  out <- list(ca0 = list(train=train.ca0, test=test.ca0),
              pco = list(train=train.pco, test=test.pco),
              cap = list(train=train.cap, test=test.cap))
  out
}

## run ranger
run_ranger_2genes <- function(dat, ntrees=500, d=d, k=2, m=2, axes=1, class="source", id="id", var_id="gene"){
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

  # find uniques
  extras <- dat.scores |> map(c("train","extra"))
  uniques <- extras |> pluck(1) |> is_unique(data=dat.test)
  
  # train model
  training <- dat.scores |> map(c("train", "training")) |> map(function(x) x |> select(-all_of(class)))
  classes <- dat.train |> pull(class)
  rf_mods <- map(training, ~ranger(classes ~ ., data=.x, oob.error = TRUE, num.trees=ntrees, respect.unordered.factors = TRUE))

  # individual tree predictions
  testing <- dat.scores |> map("test")
  tree_preds <- map2(rf_mods, testing, ~predict_by_tree(.x, .y, uniques, id=id))
  answer <- tree_preds |> map(function(x) x |> left_join(dat.test |> select(all_of(id), all_of(class))))

  list(tree_preds = answer, dat.scores=dat.scores)
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
sim_fn_2genes <- function(beta, ntrees){
  # generate data
  simdat <- gen_data_2genes(beta)
  dat = simdat$df
  d=simdat$d
  
  # run ranger predictions
  preds <- run_ranger_2genes(dat, ntrees, d, k=2, m=2, axes=1, class="source", id="id", var_id="gene")
  
  #misclassification
  mcdat <- preds$tree_preds |> bind_rows(.id = "method")
  mc <- mc_trees(mcdat, class="source")
  
  # store misclassification rates for observed alleles and unobserved alleles
  # output results
  out <- list(mc=mc, mc.extras=list(dat=dat, tree_preds=preds$tree_preds, dat.scores=preds$dat.scores))
  out
}


