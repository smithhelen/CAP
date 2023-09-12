### Simulation study CAP ###

# load libraries and functions
source("methods/libs_fns.R")
source("methods/tree_predictions.R")  # to pull out individual tree predictions
library(boot) # for inv.logit function

# start with 1 gene, 15 alleles (10 observed, 5 unobserved), 10 individuals, 2 sources
# thinking is that for observed alleles CAP will do no better than either CA or PCO
# but for unobserved alleles CAP will be better, 
# if there is a pattern in CA and strong PCO pattern in a different direction to source then CAP should be better

# set seed
set.seed(24)

## generate data
gen_data <- function(beta){
  
  # generate PCO1 and PCO2 scores
  # the direction of greatest variation is along PCO1
  PCO1 <- runif(15, -5, 5)
  PCO2 <- runif(15, -1, 1)
  
  # choose beta so when we inverse logit we get scores close to zero and close to one
  # this is the strength of differentiation along PCO2 between the sources
  # we are splitting along PCO2 whereas the direction of greatest variation is along PCO1 
  # a higher beta will have a clearer demarcation between the sources
  # this demarcation is not the same direction as the greatest spread, so won't be picked up by PCO, but should be picked up by CAP
  beta = beta # start with beta=2 # put in function as argument
  
  # generate probability scores - the probability that an allele will belong to source one.
  # a higher beta causes probs to group at zero and at one (sigmoid); a lower beta is more linear
  probs <- inv.logit(beta * PCO2)
  
  # generate counts for each allele from 10 individuals (size=10)
  # based on the probabilities generated from the position on the PCO2 axis, how many will fall in each source
  # need to do this even for the unobserved alleles so we can assign a source to the individuals later
  # this is balanced - there are 10 of each allele
  # counts of source one (s1)
  count_s1 <- map(probs, ~rbinom(n=1, size=10, p=.)) |> unlist() # find 1 random value, from a sample of 10, with probability=p
  # counts of source two (s2)
  count_s2 <- 10 - count_s1
  
  # calculate distances
  d <- data.frame(PCO1, PCO2) |> dist(diag = TRUE) |> as.matrix()
  
  # fill out individual row data ie what the original data would actually look like
  df <- data.frame(allele = factor(1:15), PCO1, PCO2, p=probs, count_s1, count_s2, observed = c(rep("Y",times=10), rep("N",times=5))) 
  dat <- df |> select(allele, PCO1, PCO2, starts_with("count"), observed) |>  
    pivot_longer(cols = starts_with("count"), names_to = "source", values_to = "n", names_prefix = "count_") |> 
    uncount(n) |> mutate(across(!where(is.double), as.factor)) |> rownames_to_column("id")
  
  out <- list(dat=dat, d=d)
  out
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
run_ranger <- function(dat, ntrees=500, d=d, k=2, m=2, axes=1, class="source", id="id", var_id="allele"){
  # prepare data - add id and split into training and test sets
  # the training data is a sample all of the individuals which have 'observed' alleles
  # the testing data has a mix of observed and unobserved alleles
  dat.train <- dat |> filter(observed == "Y") |> slice_sample(n=75)
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
sim_fn <- function(beta, ntrees){
   # generate data
   dat <- gen_data(beta)
   d <- list(dat$d)
   
   # run ranger predictions
   preds <- run_ranger(dat = dat$dat, ntrees, d=d, k=2, m=2, axes=1, class="source", id="id", var_id="allele")
   
   #misclassification
   mcdat <- preds$tree_preds |> bind_rows(.id = "method")
   mc <- mc_trees(mcdat, class="source")
   
   # store misclassification rates for observed alleles and unobserved alleles
   # output results
   out <- list(mc=mc, mc.extras=list(dat=dat, tree_preds=preds$tree_preds, dat.scores=preds$dat.scores))
   out
}



# now can change balance of alleles
# increase number of sources


