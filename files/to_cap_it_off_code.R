## R code to accompany the file "To CAP it Off - Encoding Categorical Variables using Canonical Analysis of Principal Coordinates"
# by HL Smith, PJ Biggs, NP French, ANH Smith, JC Marshall 2024

# load libraries
library(tidyverse)
library(ranger)
library(boot) # for inv.logit function
library(Manu)
library(patchwork)

# Load in our functions
source("methods/ca_unbiased.R")       # CA method with new levels scored as zero
source("methods/pco.R")               # PCO method
source("methods/cap.R")               # CAP method
source("methods/tree_predictions.R")  # to pull out individual tree predictions

# set seed
seed=456
set.seed(seed)


### Simulation study ###
## 1. generate data
# start with 1 gene, 15 levels (10 observed, 5 unobserved), 10 individuals, 2 classes
gen_data_two_class <- function(beta){
  
  # generate PCO1 and PCO2 scores
  # the direction of greatest variation is along PCO1
  # change these to change the direction to be along PCO2 for comparison
  PCO1 <- runif(15, -5, 5)
  PCO2 <- runif(15, -1, 1)
  
  # choose beta so when we inverse logit we get scores close to zero and close to one
  # this is the strength of differentiation along PCO2 between the classes
  # we are splitting along PCO2 whereas the direction of greatest variation is along PCO1 
  # a higher beta will have a clearer demarcation between the classes
  # this demarcation is not the same direction as the greatest spread, so won't be picked up by PCO, but should be picked up by CAP
  beta = beta 
  
  # generate probability scores - the probability that a level will belong to class one.
  # a higher beta causes probs to group at zero and at one (sigmoid); a lower beta is more linear
  probs <- inv.logit(beta * PCO2)
  
  # generate counts for each level from 10 individuals (size=10)
  # based on the probabilities generated from the position on the PCO2 axis, how many will fall in each class
  # need to do this even for the unobserved levels so a class can be assigned to the individuals later
  # this is balanced - there are 10 of each level
  # counts of class one (s1)
  count_s1 <- map(probs, ~rbinom(n=1, size=10, p=.)) |> unlist() # find 1 random value, from a sample of 12, with probability=p
  # counts of class two (s2)
  count_s2 <- 10 - count_s1
  
  # calculate distances
  d <- data.frame(PCO1, PCO2) |> dist(diag = TRUE) |> as.matrix()
  
  # fill out individual row data ie what the original data would actually look like
  df <- data.frame(level = factor(1:15), PCO1, PCO2, p=probs, count_s1, count_s2, observed = c(rep("Y",times=10), rep("N",times=5))) 
  dat <- df |> select(level, PCO1, PCO2, starts_with("count"), observed) |>  
    pivot_longer(cols = starts_with("count"), names_to = "class", values_to = "n", names_prefix = "count_") |> 
    uncount(n) |> rownames_to_column("id") |> mutate(across(!where(is.double), as.factor))
  
  out <- list(dat=dat, d=d, df=df)
  out
}

## 2. encode factors
# ck is the number of k axes for cap
# cm is the number of m axes for cap
# cmp is the proportion of variation captured by m axes for cap
encode <- function(Dat.train, Dat.test, d=NULL, axes=1, cmp=NULL, cm=2, ck=1, ntrees=5, id=id, class, var_id){
  train.ca0 <- prepare_training_ca0(Dat.train, starts_with(var_id), class=class, k=axes)
  test.ca0 <- prepare_test_ca0(Dat.test, train.ca0$extra, id=id)
  train.pco <- prepare_training_pco(Dat.train, starts_with(var_id), class=class, d, m=axes)
  test.pco <- prepare_test_pco(Dat.test, train.pco$extra, id=id)
  train.cap <- prepare_training_cap(Dat.train, starts_with(var_id), class=class, d, ck=ck, cm=cm, cmp=cmp, c=axes)
  test.cap <- prepare_test_cap(Dat.test, train.cap$extra, id=id)
  out <- list(ca0 = list(train=train.ca0, test=test.ca0),
              pco = list(train=train.pco, test=test.pco),
              cap = list(train=train.cap, test=test.cap))
  out
}

## 3. run ranger
run_ranger <- function(dat, ntrees=5, d=d, axes=1, cm=2, cmp=NULL, ck=1, class="class", id="id", var_id="level"){
  # prepare data - add id and split into training and test sets
  # the training data is a sample all of the individuals which have 'observed' levels
  # the testing data has a mix of observed and unobserved levels
  dat.train <- dat |> filter(observed == "Y") |> droplevels() |> slice_sample(prop = 0.75)
  dat.test <- dat |> anti_join(dat.train) |> droplevels()
  
  # set number of trees
  ntrees = ntrees
  
  # encode factors
  dat.scores <- encode(dat.train, dat.test, d=d, axes=axes, cm=cm, cmp=cmp, ck=ck, class=class, id=id, var_id=var_id, ntrees=ntrees)
  
  # find uniques
  extras <- dat.scores |> map(c("train","extra"))
  uniques <- extras |> pluck(1) |> is_unique(data=dat.test)
  
  # train model
  training <- dat.scores |> map(c("train", "training")) |> map(function(x) x |> select(-all_of(class)))
  classes <- dat.train |> pull(class)
  rf_mods <- map(training, ~ranger(classes ~ ., data=.x, oob.error = FALSE, num.trees=ntrees, sample.fraction = 1, replace = FALSE)) #, max.depth = 1
  
  # individual tree predictions
  testing <- dat.scores |> map("test")
  tree_preds <- map2(rf_mods, testing, ~predict_by_tree(.x, .y, uniques, id=id))
  answer <- tree_preds |> map(function(x) x |> left_join(dat.test |> select(all_of(id), all_of(class))))
  
  list(tree_preds = answer, dat.scores=dat.scores, dat.test=dat.test, uniques=uniques, rf_mods=rf_mods)
}

## 4. calculate misclassification rate
mc_trees <- function(dat, class){
  df <- dat |> 
    mutate(uses_unique = if_else(uses_unique == 0, "No", "Yes")) |>
    mutate(correct = !!sym(class)==prediction) |> 
    group_by(method, uses_unique, correct) |>
    summarise(n=n()) |>
    group_by(method,uses_unique) |>
    mutate(N = sum(n)) |>
    mutate(prop = n/N) |>     
    filter(correct==TRUE) |> 
    mutate(mc = 1-prop)
  df
}

# put it into a function
sim_fn <- function(beta, ntrees, axes=1, cm=2, cmp=NULL, ck=1){
  # generate data
  dat <- gen_data_two_class(beta)
  d <- list(dat$d)
  
  # run ranger predictions
  preds <- run_ranger(dat = dat$dat, ntrees, d=d, cm=cm, cmp=cmp, ck=ck, axes=axes, class="class", id="id", var_id="level")
  
  # calculate misclassification rates
  mcdat <- preds$tree_preds |> bind_rows(.id = "method")
  mc <- mc_trees(mcdat, class="class")
  
  # store misclassification rates for observed levels and unobserved levels
  # output results
  out <- list(mc=mc, mc.extras=list(dat=dat, tree_preds=preds$tree_preds, dat.scores=preds$dat.scores, dat.test=preds$dat.test, uniques=preds$uniques, rf_mods=preds$rf_mods))
  out
}

# plots for simulations
sim_plot_1a <- function(plot_dat, axis){
  plot_dat |> 
    ggplot(aes(x=.data[[axis]], y=c(rep(0, times=10)))) + 
    geom_point(aes(colour=prop_s1), size=4) +
    theme_bw(base_size = 11, base_family = 'serif') +
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.title.x = element_text(margin = margin(t=10))) + 
    scale_colour_gradientn(colours = get_pal("Pohutukawa")[c(1,4)]) +
    labs(x=paste0(axis," score"), colour="Proportion \nin class 1\n")
}

sim_plot_1b <- function(plot_dat, axis1, axis2){
  xy.limits <- range( c(pull(plot_dat, axis1), plot_dat$V2) )
  plot_dat |> 
    ggplot(aes(x=.data[[axis1]], y=V2)) + 
    geom_point(aes(colour=prop_s1), size=4) +
    theme_bw(base_size = 11, base_family = 'serif') +
    scale_colour_gradientn(colours = get_pal("Pohutukawa")[c(1,4)]) +
    scale_x_continuous(limits = xy.limits) +
    scale_y_continuous(limits = xy.limits) +
    theme(axis.title.y = element_text(margin = margin(r=5)),
          axis.title.x = element_text(margin = margin(t=10))) + 
    labs(x=paste0(axis1," score"), y=paste0(axis2," score"), colour="Proportion \nin class 1\n")
  
}

sim_plot_2 <- function(plotdat, title=NULL){
  plotdat |> mutate(method = factor(method, levels = c("ca0","pco","cap"))) |> 
    mutate(method = fct_recode(as.factor(method), `CA-unbiased`="ca0", PCO="pco", CAP="cap")) |> 
    mutate(uses_unique = factor(uses_unique, levels=c("No","Yes"))) |> 
    mutate(beta = factor(beta)) |> 
    ggplot() + 
    geom_boxplot(aes(x=method, y=mc, fill=method)) + 
    facet_grid(beta~uses_unique, 
               labeller=labeller(uses_unique=c(`No`="No absent levels \nused in prediction",
                                               `Yes`="At least one absent level \nused in prediction"),
                                 beta = c(`2`="Beta = 2",`20`="Beta = 20"))) + 
    theme_bw(base_size = 11, base_family = 'serif') + 
    labs(x="Method of encoding", y="Misclassification rate", fill="Method of encoding", title = paste0(title)) + 
    theme(strip.text.y = element_text(size = 8),
          legend.position='none',
          axis.title.x = element_text(margin = margin(t=10)), 
          axis.title.y = element_text(margin = margin(r=5))) + 
    ylim(0,1) +
    scale_fill_manual(values=manu_palettes$Kaka)
}


### Run simulations
## 1. Placement of each predictor level in PCO space vs CAP space in two dimensions (figure 2)
dat <- gen_data_two_class(beta=20)
d <- list(dat$d)
preds <- run_ranger(dat = dat$dat, ntrees=5, d=d, ck=1, axes=2, class="class", id="id", var_id="level")
Q <- preds$dat.scores$pco$train$extra$level$Q  # pco scores
pco_plot_dat <- Q |> as.data.frame() |> rename(PCO1=V1) |> rownames_to_column("level") |> 
  left_join(dat$df |> slice(1:10) |> select(level, count_s1)) |> mutate(prop_s1 = count_s1/10)
C <- preds$dat.scores$cap$train$extra$level$C_score  # cap scores
cap_plot_dat <- C |> as.data.frame() |> rename(CAP1=V1) |> rownames_to_column("level") |> 
  left_join(dat$df |> slice(1:10) |> select(level, count_s1)) |> mutate(prop_s1 = count_s1/10)

# plot results
pco1 <- sim_plot_1a(pco_plot_dat, axis="PCO1")
cap1 <- sim_plot_1a(cap_plot_dat, axis="CAP1")
pco2 <- sim_plot_1b(pco_plot_dat, "PCO1", "PCO2")
cap2 <- sim_plot_1b(cap_plot_dat, "CAP1", "CAP2")
pco2 + cap2 + pco1 + cap1 + 
  plot_layout(ncol=2, guides = "collect", heights = c(3,1)) + 
  plot_annotation(tag_levels = 'a', tag_prefix = "(", tag_suffix = ")") & 
  theme(plot.tag = element_text(size = 11, hjust = 0, vjust = 1))


## 2. Misclassification rates of trees from simulated data (figure 3)
sim_2 <- map(1:100, ~sim_fn(beta=2, ntrees=10, axes=1, cmp=100, cm=2, ck=1))
sim_20 <- map(1:100, ~sim_fn(beta=20, ntrees=10, axes=1, cmp=100, cm=2, ck=1))
mc_2 <- map(sim_2,"mc") |> bind_rows(.id="simulation") |> mutate(beta = 2)
mc_20 <- map(sim_20,"mc") |> bind_rows(.id="simulation") |> mutate(beta = 20)
cap_sim <- bind_rows(mc_2, mc_20)

# plot results
sim_plot_2(cap_sim)

