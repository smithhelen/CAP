### xgboost on SACNZ ###
# start with a few genes

# Load libraries and functions ##########################################################################################################
source("methods/libs_fns.R")
library(tidymodels)
library(finetune)
#library(xgboost)
#library(vip)

# define functions #########################################################################################################
prepare_training_pco_m <- function(data, vars, class, d, m.list, mp, residualised=NULL) {
  # pull out our var_cols and class
  var_cols <- dplyr::select(data, all_of(vars)) |> select(all_of(names(m.list)))
  classes   <- data |> pull(class)
  d.list <- d[names(m.list)]
  # iterate over the var columns and distance matrices, and convert
  prepped <- pmap(list(var_cols, d.list, m.list), factor_to_pco_score, mp) |> compact() # removes empties
  output <- map(prepped,"output")
  prepped_data <- bind_cols(data.frame(classes) |> setNames(data |> select(all_of(class)) |> colnames()),
                            map2(output, names(output), ~ .x |> set_names(paste(.y, names(.x), sep="."))))
  if(!is.null(residualised)) {
    prepped_data <- prepped_data |> bind_cols(data |> select(all_of(residualised)))
  }
  extra <- map(prepped, "extra")
  list(training = prepped_data,
       extra = extra)
}

prepare_training_cap_m <- function(data, vars, class, d, k, m.list, mp, c, residualised=NULL) {
  # pull out our var_cols and class
  var_cols <- dplyr::select(data,all_of(vars)) |> select(all_of(names(m.list)))
  classes   <- data |> pull(class)
  d.list <- d[names(m.list)]
  # iterate over the var columns and distance matrices, and convert
  prepped <- pmap(list(var_cols, d.list, m.list), factor_to_CAP_score, class = classes, k, mp, c) |> compact() # removes empties
  output <- map(prepped,"output")
  prepped_data <- bind_cols(data.frame(classes) |> setNames(data |> select(all_of(class)) |> colnames()), 
                            map2(output, names(output), ~ .x |> set_names(paste(.y, names(.x), sep="."))))
  if(!is.null(residualised)) {
    prepped_data <- prepped_data |> bind_cols(data |> select(all_of(residualised)))
  }
  extra <- map(prepped, "extra")
  list(training = prepped_data,
       extra = extra)
}

prep_data_xgb <- function(Dat.train, Dat.test, method="ca0", residualised=NULL, d=NULL, id, class, var_id, k=NULL, m.list=NULL, mp=NULL, c=NULL){
  switch(method,
         ca0 = {
           train <- prepare_training_ca0(Dat.train, starts_with(var_id), class=class, k=k, residualised=NULL)
           test <- prepare_test_ca0(Dat.test, train$extra, id=id)
         },
         pco = {
           train <- prepare_training_pco_m(Dat.train, starts_with(var_id), class=class, d=d, m.list=m.list, mp=mp, residualised=NULL)
           test <- prepare_test_pco(Dat.test, train$extra, id=id)
         },
         cap = {
           train <- prepare_training_cap_m(Dat.train, starts_with(var_id), class=class, d=d, k=k, m.list=m.list, mp=mp, c=c, residualised=NULL)
           test <- prepare_test_cap(Dat.test, train$extra, id=id)
         },
  )
  prepped <- list(train=train, test=test)
}

score_fn_m <- function(x, method, d=NULL, m.list=NULL, mp=99, k=2){
  # id the training vs test set
  df = x$data
  train.df = analysis(x)
  test.df = assessment(x)
  
  # pull out categorical column and encode
  prepped_var <- prep_data_xgb(train.df, test.df, method=method, d=d, id="LabID", class="Source", var_id="CAMP", m.list=m.list, mp=mp, k=k)
  
  # merge with numeric columns
  train_dat <- prepped_var$train$training |> bind_cols(train.df |> select(LabID))
  test_dat <- prepped_var$test |> left_join(test.df |> select(LabID, Source))
  prepped.dat <- left_join(df |> select(LabID), rbind(train_dat, test_dat)) # to get in original order (id is a factor)
  
  # return dataframe for fold and number of axes (m)
  out <- list(dat=prepped.dat, m=prepped_var$train$extra |> map("dim"))
  out
}


# Load data ##########################################################################################################
# load allele and sequence data
load("../CAP_data/data/cgMLST_dat.RData") # SACNZ cgMLST data set (cgMLST)
load("../CAP_data/data/cgMLST_seq_dat.RData") # SACNZ cgMLST sequence data set (cgMLST_seq_dat)

# distance information
load("../CAP_Data/data/list_of_distance_matrices.RData") #  distance matrices of hamming distances among unique alleles for each gene
load("../CAP_Data/data/list_of_8grams.RData")
load("Phandango/Phandango_data/list_of_distance_matrices_Phandango.RData") # hamming distances scaled to account for recombination

# list_of_distance_matrices_Phandango has some zero matrices which need to be removed
keep <- list_of_distance_matrices_Phandango[which(map(list_of_distance_matrices_Phandango, function(x) x |> sum()) != 0)] |> names()
zeros <- list_of_distance_matrices_Phandango[which(map(list_of_distance_matrices_Phandango, function(x) x |> sum()) == 0)] |> names()
list_of_distance_matrices_Phandango <- list_of_distance_matrices_Phandango[keep]

# prepare data - pull out gene names and select animal only isolates
cgMLST_dat <- cgMLST |> 
  mutate(across(everything(), factor)) |> 
  mutate(Source = factor(Source, levels = c("Beef", "Poultry", "Sheep", "Human"))) |>
  mutate(Source = fct_recode(as.factor(Source), Cattle="Beef", Chicken="Poultry")) |> 
  filter(Source != "Human") |> droplevels()

# all data
# dat <- cgMLST_dat
# dat.Ph <- cgMLST_dat |> select(!all_of(zeros)) |> droplevels()
# d <- list_of_distance_matrices
# d.8gram <- list_of_8grams
# d.Ph <- list_of_distance_matrices_Phandango

# jejuni only
details <- read.csv("../CAP_data/data/SACNZ_referencelist.csv")
dat <- cgMLST_dat |> left_join(details |> select(LabID, Species), by="LabID") |> filter(Species == "Jejuni") |> select(-c(Species, CAMP1122)) |> droplevels()
dat.Ph <- cgMLST_dat |> left_join(details |> select(LabID, Species), by="LabID") |> filter(Species == "Jejuni") |> select(-Species) |> select(!all_of(zeros)) |> droplevels()
d <- list_of_distance_matrices[colnames(dat |> select(starts_with("CAMP")))]
d.8gram <- list_of_8grams[colnames(dat |> select(starts_with("CAMP")))]
d.Ph <- list_of_distance_matrices_Phandango[colnames(dat.Ph |> select(starts_with("CAMP")))]

# subset for testing
# dat <- cgMLST_dat |> select(c(LabID, Source, sample(2:1344, size=5))) |> droplevels()
# dat.Ph <- cgMLST_dat |> select(!any_of(zeros)) |> select(c(LabID, Source, sample(2:1344, size=5))) |> droplevels()
# d <- list_of_distance_matrices[colnames(dat |> select(starts_with("CAMP")))]
# d.6gram <- list_of_6grams[colnames(dat |> select(starts_with("CAMP")))]
# d.8gram <- list_of_8grams[colnames(dat |> select(starts_with("CAMP")))]
# d.10gram <- list_of_10grams[colnames(dat |> select(starts_with("CAMP")))]
# d.12gram <- list_of_12grams[colnames(dat |> select(starts_with("CAMP")))]
# d.Ph <- list_of_distance_matrices_Phandango[colnames(dat.Ph |> select(starts_with("CAMP")))]

# seeds
runif(5,1,1000) |> round(0)

# split into folds #########################################################################################################
set.seed(483)
#nfolds=10
nfolds=5
folds <- vfold_cv(dat, strata = Source, v=nfolds)

# pull out row indices
in.id <- list()
for(i in 1:nfolds){
  in.id[[i]] <- folds$splits[[i]]$in_id
}

if(0){
  set.seed(483)
  folds.2 <- vfold_cv(dat, strata = Source, v=nfolds)
  # show they are different
  all.equal(analysis(folds$splits[[1]]), analysis(folds.2$splits[[1]]))
  for(i in 1:nfolds){
    folds.2$splits[[i]]$in_id <- in.id[[i]]
  }
  # check they are now the same
  all.equal(analysis(folds$splits[[1]]), analysis(folds.2$splits[[1]]))
}

folds.Ph <- vfold_cv(dat.Ph, strata = Source, v=nfolds)
# change indices
for(i in 1:nfolds){
  folds.Ph$splits[[i]]$in_id <- in.id[[i]]
}

# score each fold #########################################################################################################
ca0_folds <- folds
for(i in 1:nfolds){
  x <- ca0_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="ca0", d=NULL)
  ca0_folds$splits[[i]]$data <- prepped.dat$dat
}

#find common variables
cols <- vector(mode="list", length = nfolds)
for(i in 1:nfolds){
  cols[[i]] <- ca0_folds$splits[[i]]$data |> colnames()
}
cols <- reduce(cols, intersect)
for(i in 1:nfolds){
  ca0_folds$splits[[i]]$data <- ca0_folds$splits[[i]]$data |> select(all_of(cols))
}


pco_folds <- folds
# find the maximum number of pco axes for each gene

# define list of m=null so that the number of axes will be determined by propG/mp
m.null <- list(NULL) |> rep(length(d))
names(m.null) <- names(d)

# create empty list for mp=95%
m <- list()
for(i in 1:nfolds){
  x <- pco_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d, m.null, mp=95)$m
}
m.95 <- map(m, unlist) |> bind_rows() |> map(max, na.rm=TRUE)

# create empty list for mp=99%
m <- list()
for(i in 1:nfolds){
  x <- pco_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d, m.null, mp=99)$m
}
m.99 <- map(m, unlist) |> bind_rows() |> map(max, na.rm=TRUE)

# now use m (number of axes) rather than mp (prop var)
for(i in 1:nfolds){
  x <- pco_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="pco", d=d, m.list=m.95, mp=NULL)
  pco_folds$splits[[i]]$data <- prepped.dat$dat
}
#pco_folds$splits |> map(function(x) x$data |> ncol()) 
#pco_folds$splits |> map(function(x) x$data |> colnames()) 

#find common variables
cols <- vector(mode="list", length = nfolds)
for(i in 1:nfolds){
  cols[[i]] <- pco_folds$splits[[i]]$data |> colnames()
}
cols <- reduce(cols, intersect)
for(i in 1:nfolds){
  pco_folds$splits[[i]]$data <- pco_folds$splits[[i]]$data |> select(all_of(cols))
  }

# cap_folds <- folds # use m.99 as m.95 has a lot of single dimension variables
# for(i in 1:nfolds){
#   x <- cap_folds$splits[[i]]
#   prepped.dat <- score_fn_m(x, method="cap", m.list=m.99, mp=NULL, d=d)
#   cap_folds$splits[[i]]$data <- prepped.dat$dat
# }
cap_folds <- folds # m.95 use for jejuni
for(i in 1:nfolds){
  x <- cap_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="cap", m.list=m.95, mp=NULL, d=d)
  cap_folds$splits[[i]]$data <- prepped.dat$dat
}
#cap_folds$splits |> map(function(x) x$data |> ncol()) 
#cap_folds$splits |> map(function(x) x$data |> colnames()) 
#find common variables
cols <- vector(mode="list", length = nfolds)
for(i in 1:nfolds){
  cols[[i]] <- cap_folds$splits[[i]]$data |> colnames()
}
cols <- reduce(cols, intersect)
for(i in 1:nfolds){
  cap_folds$splits[[i]]$data <- cap_folds$splits[[i]]$data |> select(all_of(cols))
}

pco_8gram_folds <- folds
# find the maximum number of pco axes for each gene

# create empty list for mp=95%
m <- list()
for(i in 1:nfolds){
  x <- pco_8gram_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d.8gram, m.null, mp=95)$m
}
m.95 <- map(m, unlist) |> bind_rows() |> map(max, na.rm=TRUE)
# create empty list for mp=99%
m <- list()
for(i in 1:nfolds){
  x <- pco_8gram_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d.8gram, m.null, mp=99)$m
}
m.99 <- map(m, unlist) |> bind_rows() |> map(max, na.rm=TRUE)

# now use m (number of axes) rather than mp (prop var)
for(i in 1:nfolds){
  x <- pco_8gram_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="pco", d=d.8gram, m.list=m.95, mp=NULL)
  pco_8gram_folds$splits[[i]]$data <- prepped.dat$dat
}
#pco_8gram_folds$splits |> map(function(x) x$data |> ncol()) 
#pco_8gram_folds$splits |> map(function(x) x$data |> colnames()) 
cols <- vector(mode="list", length = nfolds)
for(i in 1:nfolds){
  cols[[i]] <- pco_8gram_folds$splits[[i]]$data |> colnames()
}
cols <- reduce(cols, intersect)
for(i in 1:nfolds){
    pco_8gram_folds$splits[[i]]$data <- pco_8gram_folds$splits[[i]]$data |> select(all_of(cols))
}



# cap_8gram_folds <- folds # use m.99 as m.95 has a lot of single dimension variables
# for(i in 1:nfolds){
#   x <- cap_8gram_folds$splits[[i]]
#   prepped.dat <- score_fn_m(x, method="cap", m.list=m.99, mp=NULL, d=d.8gram, k=2)
#   cap_8gram_folds$splits[[i]]$data <- prepped.dat$dat
# }
cap_8gram_folds <- folds # m.95 use for jejuni
for(i in 1:nfolds){
  x <- cap_8gram_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="cap", m.list=m.95, mp=NULL, d=d.8gram)
  cap_8gram_folds$splits[[i]]$data <- prepped.dat$dat
}
#cap_8gram_folds$splits |> map(function(x) x$data |> ncol()) 
#cap_8gram_folds$splits |> map(function(x) x$data |> ncol()) |> which.min()
#cap_8gram_folds$splits |> map(function(x) x$data |> colnames()) 
cols <- vector(mode="list", length = nfolds)
for(i in 1:nfolds){
  cols[[i]] <- cap_8gram_folds$splits[[i]]$data |> colnames()
}
cols <- reduce(cols, intersect)
for(i in 1:nfolds){
  cap_8gram_folds$splits[[i]]$data <- cap_8gram_folds$splits[[i]]$data |> select(all_of(cols))
}



pco_Phand_folds <- folds.Ph
# find the maximum number of pco axes for each gene
m.null.Ph <- list(NULL) |> rep(length(d.Ph))
names(m.null.Ph) <- names(d.Ph)

# create empty list for mp=95%
m <- list()
for(i in 1:nfolds){
  x <- pco_Phand_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d.Ph, m.null.Ph, mp=95)$m
}
m.95 <- map(m, unlist) |> bind_rows() |> map(max, na.rm=TRUE)
length(m.95)
# create empty list for mp=99%
m <- list()
for(i in 1:nfolds){
  x <- pco_Phand_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d.Ph, m.null.Ph, mp=99)$m
}
m.99 <- map(m, unlist) |> bind_rows() |> map(max, na.rm=TRUE)

# now use m (number of axes) rather than mp (prop var)
for(i in 1:nfolds){
  x <- pco_Phand_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="pco", d=d.Ph, m.list=m.95, mp=NULL)
  pco_Phand_folds$splits[[i]]$data <- prepped.dat$dat
}
#pco_Phand_folds$splits |> map(function(x) x$data |> ncol()) 
#pco_Phand_folds$splits |> map(function(x) x$data |> colnames()) 
cols <- vector(mode="list", length = nfolds)
for(i in 1:nfolds){
  cols[[i]] <- pco_Phand_folds$splits[[i]]$data |> colnames()
}
cols <- reduce(cols, intersect)
for(i in 1:nfolds){
    pco_Phand_folds$splits[[i]]$data <- pco_Phand_folds$splits[[i]]$data |> select(all_of(cols))
  }


# cap_Phand_folds <- folds.Ph # use m.99 as m.95 has a lot of single dimension variables
# for(i in 1:nfolds){
#   x <- cap_Phand_folds$splits[[i]]
#   prepped.dat <- score_fn_m(x, method="cap", m.list=m.99, mp=NULL, d=d.Ph, k=2)
#   cap_Phand_folds$splits[[i]]$data <- prepped.dat$dat
# }
cap_Phand_folds <- folds.Ph # m.95 use for jejuni
for(i in 1:nfolds){
  x <- cap_Phand_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="cap", m.list=m.95, mp=NULL, d=d.Ph)
  cap_Phand_folds$splits[[i]]$data <- prepped.dat$dat
}
#cap_Phand_folds$splits |> map(function(x) x$data |> ncol()) 
#cap_Phand_folds$splits |> map(function(x) x$data |> ncol()) |> which.min()
#cap_Phand_folds$splits |> map(function(x) x$data |> colnames()) 
cols <- vector(mode="list", length = nfolds)
for(i in 1:nfolds){
  cols[[i]] <- cap_Phand_folds$splits[[i]]$data |> colnames()
}
cols <- reduce(cols, intersect)
for(i in 1:nfolds){
    cap_Phand_folds$splits[[i]]$data <- cap_Phand_folds$splits[[i]]$data |> select(all_of(cols))
}





# define xgboost model, recipe and workflow #########################################################################################################
xg_spec <- 
  boost_tree(trees = tune(), #default 15
             min_n = tune(), #default 1
             mtry = tune(), #default 
             #tree_depth = tune(), #default 6           
             #learn_rate = 0.01
             learn_rate = tune()) |>
  set_engine("xgboost") |> 
  set_mode("classification")

xg_recipe.ca0 <- recipe(Source ~ ., data = ca0_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.ca0 <- workflow(xg_recipe.ca0, xg_spec)

xg_recipe.pco <- recipe(Source ~ ., data = pco_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.pco <- workflow(xg_recipe.pco, xg_spec)

xg_recipe.cap <- recipe(Source ~ ., data = cap_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.cap <- workflow(xg_recipe.cap, xg_spec)

xg_recipe.pco.8g <- recipe(Source ~ ., data = pco_8gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.pco.8g <- workflow(xg_recipe.pco.8g, xg_spec)

xg_recipe.cap.8g <- recipe(Source ~ ., data = cap_8gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.cap.8g <- workflow(xg_recipe.cap.8g, xg_spec)

xg_recipe.pco.Ph <- recipe(Source ~ ., data = pco_Phand_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.pco.Ph <- workflow(xg_recipe.pco.Ph, xg_spec)

xg_recipe.cap.Ph <- recipe(Source ~ ., data = cap_Phand_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.cap.Ph <- workflow(xg_recipe.cap.Ph, xg_spec)


# tune xgboost parameters ##############################################################################################################
doParallel::registerDoParallel()
set.seed(568)
xg_resamples.ca0 <- tune_race_anova(
  xg_workflow.ca0,
  resamples = ca0_folds,
  grid = 15,
  control = control_race(verbose_elim = TRUE)
)

xg_resamples.pco <- tune_race_anova(
  xg_workflow.pco,
  resamples = pco_folds,
  grid = 15,
  control = control_race(verbose_elim = TRUE)
)

xg_resamples.cap <- tune_race_anova(
  xg_workflow.cap,
  resamples = cap_folds,
  grid = 15,
  control = control_race(verbose_elim = TRUE)
)

xg_resamples.pco.8g <- tune_race_anova(
  xg_workflow.pco.8g,
  resamples = pco_8gram_folds,
  grid = 15,
  control = control_race(verbose_elim = TRUE)
)

xg_resamples.cap.8g <- tune_race_anova(
  xg_workflow.cap.8g,
  resamples = cap_8gram_folds,
  grid = 15,
  control = control_race(verbose_elim = TRUE)
)

xg_resamples.pco.Ph <- tune_race_anova(
  xg_workflow.pco.Ph,
  resamples = pco_Phand_folds,
  grid = 15,
  control = control_race(verbose_elim = TRUE)
)

xg_resamples.cap.Ph <- tune_race_anova(
  xg_workflow.cap.Ph,
  resamples = cap_Phand_folds,
  grid = 15,
  control = control_race(verbose_elim = TRUE)
)



# find best parameters ##############################################################################################################
#collect_metrics(xg_resamples.ca0)
#plot_race(xg_resamples.ca0)
#show_best(xg_resamples.ca0, "accuracy")

best_ca0 <- select_best(xg_resamples.ca0, "accuracy")
best_pco <- select_best(xg_resamples.pco, "accuracy")
best_cap <- select_best(xg_resamples.cap, "accuracy")
best_pco.8g <- select_best(xg_resamples.pco.8g, "accuracy")
best_cap.8g <- select_best(xg_resamples.cap.8g, "accuracy")
best_pco.Ph <- select_best(xg_resamples.pco.Ph, "accuracy")
best_cap.Ph <- select_best(xg_resamples.cap.Ph, "accuracy")


# split into folds ##############################################################################################################
set.seed(458)
folds <- vfold_cv(dat, strata = Source, v=nfolds)

# pull out row indices
in.id <- list()
for(i in 1:nfolds){
  in.id[[i]] <- folds$splits[[i]]$in_id
}

folds.Ph <- vfold_cv(dat.Ph, strata = Source, v=nfolds)
# change indices
for(i in 1:nfolds){
  folds.Ph$splits[[i]]$in_id <- in.id[[i]]
}

# score each fold ##############################################################################################################
ca0_folds <- folds
for(i in 1:nfolds){
  x <- ca0_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="ca0", d=NULL)
  ca0_folds$splits[[i]]$data <- prepped.dat$dat
}

cols <- vector(mode="list", length = nfolds)
for(i in 1:nfolds){
  cols[[i]] <- ca0_folds$splits[[i]]$data |> colnames()
}
cols <- reduce(cols, intersect)
for(i in 1:nfolds){
    ca0_folds$splits[[i]]$data <- ca0_folds$splits[[i]]$data |> select(all_of(cols))
  }



pco_folds <- folds
# find the maximum number of pco axes for each gene

# define list of m=null so that the number of axes will be determined by propG/mp
m.null <- list(NULL) |> rep(length(d))
names(m.null) <- names(d)

# create empty list for mp=95%
m <- list()
for(i in 1:nfolds){
  x <- pco_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d, m.null, mp=95)$m
}
m.95 <- map(m, unlist) |> bind_rows() |> map(max, na.rm=TRUE)

# create empty list for mp=99%
m <- list()
for(i in 1:nfolds){
  x <- pco_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d, m.null, mp=99)$m
}
m.99 <- map(m, unlist) |> bind_rows() |> map(max, na.rm=TRUE)

# now use m (number of axes) rather than mp (prop var)
for(i in 1:nfolds){
  x <- pco_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="pco", d=d, m.list=m.95, mp=NULL)
  pco_folds$splits[[i]]$data <- prepped.dat$dat
}
#pco_folds$splits |> map(function(x) x$data |> ncol()) 
#pco_folds$splits |> map(function(x) x$data |> colnames()) 
cols <- vector(mode="list", length = nfolds)
for(i in 1:nfolds){
  cols[[i]] <- pco_folds$splits[[i]]$data |> colnames()
}
cols <- reduce(cols, intersect)
for(i in 1:nfolds){
    pco_folds$splits[[i]]$data <- pco_folds$splits[[i]]$data |> select(all_of(cols))
  }



# cap_folds <- folds # use m.99 as m.95 has a lot of single dimension variables
# for(i in 1:nfolds){
#   x <- cap_folds$splits[[i]]
#   prepped.dat <- score_fn_m(x, method="cap", m.list=m.99, mp=NULL, d=d)
#   cap_folds$splits[[i]]$data <- prepped.dat$dat
# }
cap_folds <- folds # m.95 for jejuni
for(i in 1:nfolds){
  x <- cap_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="cap", m.list=m.95, mp=NULL, d=d)
  cap_folds$splits[[i]]$data <- prepped.dat$dat
}
#cap_folds$splits |> map(function(x) x$data |> ncol()) 
#cap_folds$splits |> map(function(x) x$data |> colnames()) 
cols <- vector(mode="list", length = nfolds)
for(i in 1:nfolds){
  cols[[i]] <- cap_folds$splits[[i]]$data |> colnames()
}
cols <- reduce(cols, intersect)
for(i in 1:nfolds){
    cap_folds$splits[[i]]$data <- cap_folds$splits[[i]]$data |> select(all_of(cols))
  }



pco_8gram_folds <- folds
# find the maximum number of pco axes for each gene

# create empty list for mp=95%
m <- list()
for(i in 1:nfolds){
  x <- pco_8gram_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d.8gram, m.null, mp=95)$m
}
m.95 <- map(m, unlist) |> bind_rows() |> map(max, na.rm=TRUE)
# create empty list for mp=99%
m <- list()
for(i in 1:nfolds){
  x <- pco_8gram_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d.8gram, m.null, mp=99)$m
}
m.99 <- map(m, unlist) |> bind_rows() |> map(max, na.rm=TRUE)

# now use m (number of axes) rather than mp (prop var)
for(i in 1:nfolds){
  x <- pco_8gram_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="pco", d=d.8gram, m.list=m.95, mp=NULL)
  pco_8gram_folds$splits[[i]]$data <- prepped.dat$dat
}
#pco_8gram_folds$splits |> map(function(x) x$data |> ncol()) 
#pco_8gram_folds$splits |> map(function(x) x$data |> colnames()) 
cols <- vector(mode="list", length = nfolds)
for(i in 1:nfolds){
  cols[[i]] <- pco_8gram_folds$splits[[i]]$data |> colnames()
}
cols <- reduce(cols, intersect)
for(i in 1:nfolds){
    pco_8gram_folds$splits[[i]]$data <- pco_8gram_folds$splits[[i]]$data |> select(all_of(cols))
  }


# cap_8gram_folds <- folds # use m.99 as m.95 has a lot of single dimension variables
# for(i in 1:nfolds){
#   x <- cap_8gram_folds$splits[[i]]
#   prepped.dat <- score_fn_m(x, method="cap", m.list=m.99, mp=NULL, d=d.8gram)
#   cap_8gram_folds$splits[[i]]$data <- prepped.dat$dat
# }
cap_8gram_folds <- folds # m.95 for jejuni
for(i in 1:nfolds){
  x <- cap_8gram_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="cap", m.list=m.95, mp=NULL, d=d.8gram)
  cap_8gram_folds$splits[[i]]$data <- prepped.dat$dat
}
#cap_8gram_folds$splits |> map(function(x) x$data |> ncol()) 
#cap_8gram_folds$splits |> map(function(x) x$data |> colnames()) 
cols <- vector(mode="list", length = nfolds)
for(i in 1:nfolds){
  cols[[i]] <- cap_8gram_folds$splits[[i]]$data |> colnames()
}
cols <- reduce(cols, intersect)
cap_8gram_folds$splits |> map(function(x) x$data |> ncol()) 
  for(i in 1:nfolds){
    cap_8gram_folds$splits[[i]]$data <- cap_8gram_folds$splits[[i]]$data |> select(all_of(cols))
  }




pco_Phand_folds <- folds.Ph
# find the maximum number of pco axes for each gene

# create empty list for mp=95%
m <- list()
for(i in 1:nfolds){
  x <- pco_Phand_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d.Ph, m.null.Ph, mp=95)$m
}
m.95 <- map(m, unlist) |> bind_rows() |> map(max, na.rm=TRUE)
# create empty list for mp=99%
m <- list()
for(i in 1:nfolds){
  x <- pco_Phand_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d.Ph, m.null.Ph, mp=99)$m
}
m.99 <- map(m, unlist) |> bind_rows() |> map(max, na.rm=TRUE)

# now use m (number of axes) rather than mp (prop var)
for(i in 1:nfolds){
  x <- pco_Phand_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="pco", d=d.Ph, m.list=m.95, mp=NULL)
  pco_Phand_folds$splits[[i]]$data <- prepped.dat$dat
}
#pco_Phand_folds$splits |> map(function(x) x$data |> ncol()) 
#pco_Phand_folds$splits |> map(function(x) x$data |> colnames()) 
cols <- vector(mode="list", length = nfolds)
for(i in 1:nfolds){
  cols[[i]] <- pco_Phand_folds$splits[[i]]$data |> colnames()
}
cols <- reduce(cols, intersect)
for(i in 1:nfolds){
    pco_Phand_folds$splits[[i]]$data <- pco_Phand_folds$splits[[i]]$data |> select(all_of(cols))
  }


# 
# cap_Phand_folds <- folds.Ph # use m.99 as m.95 has a lot of single dimension variables
# for(i in 1:nfolds){
#   x <- cap_Phand_folds$splits[[i]]
#   prepped.dat <- score_fn_m(x, method="cap", m.list=m.99, mp=NULL, d=d.Ph)
#   cap_Phand_folds$splits[[i]]$data <- prepped.dat$dat
# }
cap_Phand_folds <- folds.Ph # m.95 for jejuni
for(i in 1:nfolds){
  x <- cap_Phand_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="cap", m.list=m.95, mp=NULL, d=d.Ph)
  cap_Phand_folds$splits[[i]]$data <- prepped.dat$dat
}
#cap_Phand_folds$splits |> map(function(x) x$data |> ncol()) 
#cap_Phand_folds$splits |> map(function(x) x$data |> colnames()) 
cols <- vector(mode="list", length = nfolds)
for(i in 1:nfolds){
  cols[[i]] <- cap_Phand_folds$splits[[i]]$data |> colnames()
}
cols <- reduce(cols, intersect)
for(i in 1:nfolds){
    cap_Phand_folds$splits[[i]]$data <- cap_Phand_folds$splits[[i]]$data |> select(all_of(cols))
  }



# update model specs ############################################################################################################
# need to update the recipe as well as the columns might be different in the new data folds
xg_recipe.ca0 <- recipe(Source ~ ., data = ca0_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.ca0 <- workflow(xg_recipe.ca0, xg_spec)
xg_final_workflow.ca0 <- finalize_workflow(xg_workflow.ca0, best_ca0)

xg_recipe.pco <- recipe(Source ~ ., data = pco_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.pco <- workflow(xg_recipe.pco, xg_spec)
xg_final_workflow.pco <- finalize_workflow(xg_workflow.pco, best_pco)

xg_recipe.cap <- recipe(Source ~ ., data = cap_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.cap <- workflow(xg_recipe.cap, xg_spec)
xg_final_workflow.cap <- finalize_workflow(xg_workflow.cap, best_cap)

xg_recipe.pco.8g <- recipe(Source ~ ., data = pco_8gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.pco.8g <- workflow(xg_recipe.pco.8g, xg_spec)
xg_final_workflow.pco.8g <- finalize_workflow(xg_workflow.pco.8g, best_pco.8g)

xg_recipe.cap.8g <- recipe(Source ~ ., data = cap_8gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.cap.8g <- workflow(xg_recipe.cap.8g, xg_spec)
xg_final_workflow.cap.8g <- finalize_workflow(xg_workflow.cap.8g, best_cap.8g)

xg_recipe.pco.Ph <- recipe(Source ~ ., data = pco_Phand_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.pco.Ph <- workflow(xg_recipe.pco.Ph, xg_spec)
xg_final_workflow.pco.Ph <- finalize_workflow(xg_workflow.pco.Ph, best_pco.Ph)

xg_recipe.cap.Ph <- recipe(Source ~ ., data = cap_Phand_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.cap.Ph <- workflow(xg_recipe.cap.Ph, xg_spec)
xg_final_workflow.cap.Ph <- finalize_workflow(xg_workflow.cap.Ph, best_cap.Ph)


xg_workflows <- list(ca0=xg_final_workflow.ca0,
                     pco=xg_final_workflow.pco,
                     cap=xg_final_workflow.cap,
                     pco_8g=xg_final_workflow.pco.8g,
                     cap_8g=xg_final_workflow.cap.8g,
                     pco_Ph=xg_final_workflow.pco.Ph,
                     cap_Ph=xg_final_workflow.cap.Ph)
#save(xg_workflows, file = "xgb_workflows_SACNZ.Rdata")

# run xgboost resamples ##############################################################################################################
set.seed(659)
xg_rs.ca0 <- fit_resamples(
  xg_final_workflow.ca0,
  resamples = ca0_folds,
  control = control_resamples(save_pred = TRUE)
)
#collect_metrics(xg_rs.ca0)

xg_rs.pco <- fit_resamples(
  xg_final_workflow.pco,
  resamples = pco_folds,
  control = control_resamples(save_pred = TRUE)
)
#collect_metrics(xg_rs.pco)

xg_rs.cap <- fit_resamples(
  xg_final_workflow.cap,
  resamples = cap_folds,
  control = control_resamples(save_pred = TRUE)
)
#collect_metrics(xg_rs.cap)

xg_rs.pco.8g <- fit_resamples(
  xg_final_workflow.pco.8g,
  resamples = pco_8gram_folds,
  control = control_resamples(save_pred = TRUE)
)
#collect_metrics(xg_rs.pco.8g)

xg_rs.cap.8g <- fit_resamples(
  xg_final_workflow.cap.8g,
  resamples = cap_8gram_folds,
  control = control_resamples(save_pred = TRUE)
)
#collect_metrics(xg_rs.cap.8g)

xg_rs.pco.Ph <- fit_resamples(
  xg_final_workflow.pco.Ph,
  resamples = pco_Phand_folds,
  control = control_resamples(save_pred = TRUE)
)
#collect_metrics(xg_rs.pco.Ph)

xg_rs.cap.Ph <- fit_resamples(
  xg_final_workflow.cap.Ph,
  resamples = cap_Phand_folds,
  control = control_resamples(save_pred = TRUE)
)
#collect_metrics(xg_rs.cap.Ph)


xgb_metrics <- list(ca0=xg_rs.ca0, 
                    pco=xg_rs.pco, 
                    cap=xg_rs.cap, 
                    pco_8g=xg_rs.pco.8g, 
                    cap_8g=xg_rs.cap.8g, 
                    pco_Ph=xg_rs.pco.Ph, 
                    cap_Ph=xg_rs.cap.Ph)
#save(xgb_metrics, file = "xgb_metrics_SACNZ.Rdata")

# define rf model, recipe and workflow ##############################################################################################################
rf_spec <- 
  rand_forest(trees = 500)  |>  
  set_engine("ranger", respect.unordered.factors = TRUE) |> 
  set_mode("classification")

rf_recipe.ca0 <- recipe(Source ~ ., data = ca0_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
rf_workflow.ca0 <- workflow(rf_recipe.ca0, rf_spec)

rf_recipe.pco <- recipe(Source ~ ., data = pco_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
rf_workflow.pco <- workflow(rf_recipe.pco, rf_spec)

rf_recipe.cap <- recipe(Source ~ ., data = cap_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
rf_workflow.cap <- workflow(rf_recipe.cap, rf_spec)

rf_recipe.pco.8g <- recipe(Source ~ ., data = pco_8gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
rf_workflow.pco.8g <- workflow(rf_recipe.pco.8g, rf_spec)

rf_recipe.cap.8g <- recipe(Source ~ ., data = cap_8gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
rf_workflow.cap.8g <- workflow(rf_recipe.cap.8g, rf_spec)

rf_recipe.pco.Ph <- recipe(Source ~ ., data = pco_Phand_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
rf_workflow.pco.Ph <- workflow(rf_recipe.pco.Ph, rf_spec)

rf_recipe.cap.Ph <- recipe(Source ~ ., data = cap_Phand_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
rf_workflow.cap.Ph <- workflow(rf_recipe.cap.Ph, rf_spec)

# run rf resamples ##############################################################################################################

set.seed(847)
rf_rs.ca0 <- fit_resamples(
  rf_workflow.ca0,
  resamples = ca0_folds,
  control = control_resamples(save_pred = TRUE)
)
#collect_metrics(rf_rs.ca0)

rf_rs.pco <- fit_resamples(
  rf_workflow.pco,
  resamples = pco_folds,
  control = control_resamples(save_pred = TRUE)
)
#collect_metrics(rf_rs.pco)

rf_rs.cap <- fit_resamples(
  rf_workflow.cap,
  resamples = cap_folds,
  control = control_resamples(save_pred = TRUE)
)
#collect_metrics(rf_rs.cap)

rf_rs.pco.8g <- fit_resamples(
  rf_workflow.pco.8g,
  resamples = pco_8gram_folds,
  control = control_resamples(save_pred = TRUE)
)
#collect_metrics(rf_rs.pco.8g)

rf_rs.cap.8g <- fit_resamples(
  rf_workflow.cap.8g,
  resamples = cap_8gram_folds,
  control = control_resamples(save_pred = TRUE)
)
#collect_metrics(rf_rs.cap.8g)

rf_rs.pco.Ph <- fit_resamples(
  rf_workflow.pco.Ph,
  resamples = pco_Phand_folds,
  control = control_resamples(save_pred = TRUE)
)
#collect_metrics(rf_rs.pco.Ph)

rf_rs.cap.Ph <- fit_resamples(
  rf_workflow.cap.Ph,
  resamples = cap_Phand_folds,
  control = control_resamples(save_pred = TRUE)
)
#collect_metrics(rf_rs.cap.Ph)

rf_metrics <- list(ca0=rf_rs.ca0,
                   pco=rf_rs.pco,
                   cap=rf_rs.cap,
                   pco_8g=rf_rs.pco.8g,
                   cap_8g=rf_rs.cap.8g,
                   pco_Ph=rf_rs.pco.Ph,
                   cap_Ph=rf_rs.cap.Ph
                   )
#save(rf_metrics, file = "test_rf_metrics_SACNZ.Rdata")

all_metrics <- list(ca0=list(rf=rf_rs.ca0, xg=xg_rs.ca0),
                    pco=list(rf=rf_rs.pco, xg=xg_rs.pco),
                    cap=list(rf=rf_rs.cap,xg=xg_rs.cap),
                    pco_8g=list(rf=rf_rs.pco.8g, xg=xg_rs.pco.8g),
                    cap_8g=list(rf=rf_rs.cap.8g,xg=xg_rs.cap.8g),
                    pco_Ph=list(rf=rf_rs.pco.Ph, xg=xg_rs.pco.Ph),
                    cap_Ph=list(rf=rf_rs.cap.Ph,xg=xg_rs.cap.Ph)
                    )
save(all_metrics, file = "all_metrics_SACNZ_5folds_jejuni_v2.Rdata")
#load("all_metrics_SACNZ_5folds_v1.Rdata")
#load("all_metrics_SACNZ_5folds_v2.Rdata")
#load("all_metrics_SACNZ_v1.Rdata")
#load("all_metrics_SACNZ_v2.Rdata")
#load("all_metrics_SACNZ_v3.Rdata")


# plot ##############################################################################################################
nfolds = length(all_metrics$cap$rf$splits)

p <- all_metrics |> map_df(function(x) x |> map_df(collect_metrics, .id = "model"), .id = "method" ) |> 
  filter(.metric == "accuracy") |> 
  mutate(dist = case_when(method |> str_detect("8g") ~ 15,
                          method |> str_detect("ca0") ~ 1,
                          method |> str_detect("Ph") ~ 17,
                          TRUE ~ 19),
         method2 = factor(method),
         method = factor(str_sub(method2, 1,3), ordered=TRUE, levels=c("ca0", "pco", "cap"))) |> 
  select(method2, method, model, mean, std_err, dist)

p |> 
  ggplot(aes(x=method, y=mean)) +
  geom_errorbar(aes(col=method, group=as.factor(dist), ymin=mean-std_err, ymax=mean+std_err), width=0.2, position = position_dodge(0.5)) +
  geom_point(aes(col=method, shape=dist, group=as.factor(dist)), size=3, stroke=1.5, position = position_dodge(0.5)) +
  #ylim(0.55,0.9) +
  scale_shape_identity(guide = "legend", breaks=c(1, 15,17,19), labels=c("None","8grams","Phandango","Hamming")) +
  scale_colour_manual(values = c("#115896", "#4C4C53", "#BA2F00")) +
  guides(col = guide_legend("Method of encoding"),
         shape = guide_legend("Distance measure")) +
  labs(x = "Method of encoding", y = "Classification success", title = paste0("xgboost and rf results, SACNZ sample data, ", nfolds ," cv folds")) +
  theme_bw() +
  facet_wrap(~model, labeller = labeller(model=c("rf" = "Random Forest", "xg" = "xgboost")))

# by fold
p2 <- bind_rows(
  rf = all_metrics |> map_df(function(x) x$rf |> pull(.metrics) |> map_df(function(x) x |> filter(.metric=="accuracy") |> select(.estimate), .id="fold"), .id = "method" ),
  xg = all_metrics |> map_df(function(x) x$xg |> pull(.metrics) |> map_df(function(x) x |> filter(.metric=="accuracy") |> select(.estimate), .id="fold"), .id = "method" ), .id="model") |> 
  mutate(dist = case_when(method |> str_detect("8g") ~ 15,
                          method |> str_detect("ca0") ~ 1,
                          method |> str_detect("Ph") ~ 17,
                          TRUE ~ 19),
         method2 = factor(method),
         method = factor(str_sub(method2, 1,3), ordered=TRUE, levels=c("ca0", "pco", "cap")))

p2 |> 
  ggplot(aes(x=method, y=.estimate)) +
  geom_point(aes(col=method, shape=dist, group=as.factor(dist)), size=3, stroke=1.5, position = position_dodge(0.5)) +
  scale_shape_identity(guide = "legend", breaks=c(1, 15,17,19), labels=c("None","8grams","Phandango","Hamming")) +
  stat_summary(aes(group=as.factor(dist)), col="black", geom = "point", fun = "mean", size = 4, stroke=1.5, shape = 4, position = position_dodge(0.5), show.legend = FALSE) +
  ylim(0.,1.0) +
  scale_colour_manual(values = c("#086abf", "#6c7a87", "#da400b")) +
  guides(col = guide_legend("Method of encoding"),
         shape = guide_legend("Distance measure")) +
  labs(x = "Method of encoding", y = "Classification success", title = paste0("xgboost and rf results, SACNZ sample data, ", nfolds, " cv folds")) +
  theme_bw() +
  facet_wrap(~model, labeller = labeller(model=c("rf" = "Random Forest", "xg" = "xgboost")))

p2 |> 
  ggplot(aes(x=method, y=.estimate)) +
  geom_text(aes(label=fold, colour=as.factor(dist), group=as.factor(dist)), size=4, position = position_dodge(0.5))  + 
  geom_point(aes(colour=as.factor(dist)), alpha = 0, size=3) + 
  stat_summary(aes(group=as.factor(dist)), col="black", geom = "point", fun = "mean", size = 4, stroke=1.5, shape = 4, position = position_dodge(0.5), show.legend = FALSE) +
 # ylim(0.50,0.90) +
  scale_colour_manual(values = c("#086abf", "#6c7a87", "#da400b", "#3da572"), guide = "legend", labels=c("None","8grams","Phandango","Hamming")) +
  guides(col = guide_legend(title="Distance measure", override.aes=aes(label="", alpha=1))) +
  labs(x = "Method of encoding", y = "Classification success", title = paste0("xgboost and rf results, SACNZ sample data, ", nfolds, " cv folds")) +
  theme_bw() +
  facet_wrap(~model, labeller = labeller(model=c("rf" = "Random Forest", "xg" = "xgboost")))


# by class
p3 <- bind_rows(
  rf = all_metrics |> map_df(function(x) x$rf$.predictions |> reduce(rbind) |> select(.pred_class, Source) |> group_by(Source) |> 
                               summarise(N=n(), n=sum(.pred_class==Source), p=n/N), .id="method"),
  xg = all_metrics |> map_df(function(x) x$xg$.predictions |> reduce(rbind) |> select(.pred_class, Source) |> group_by(Source) |> 
                               summarise(N=n(), n=sum(.pred_class==Source), p=n/N), .id="method"), .id = "model"
  ) |> 
  mutate(dist = case_when(method |> str_detect("8g") ~ 15,
                          method |> str_detect("ca0") ~ 1,
                          method |> str_detect("Ph") ~ 17,
                          TRUE ~ 19),
         method2 = factor(method),
         method = factor(str_sub(method2, 1,3), ordered=TRUE, levels=c("ca0", "pco", "cap")))
  
p3 |> 
  ggplot(aes(x=method, y=p)) +
  geom_point(aes(col=Source, shape=dist, group=as.factor(dist)), size=3, stroke=1.5, position = position_dodge(0.5)) +
  scale_shape_identity(guide = "legend", breaks=c(15,17,19), labels=c("8grams","Phandango","Hamming")) +
  #stat_summary(aes(group=as.factor(dist), col=method), geom = "point", fun = "mean", size = 5, shape = 4, position = position_dodge(0.5)) +
  ylim(0,1) +
  scale_colour_manual(values = c("#115896", "#4C4C53", "#BA2F00")) +
  guides(col = guide_legend("Source")) +
  labs(x = "Method of encoding", y = "Classification success", title = paste0("xgboost and rf results, SACNZ sample data, ", nfolds, " cv folds")) +
  theme_bw() +
  facet_wrap(~model, labeller = labeller(model=c("rf" = "Random Forest", "xg" = "xgboost")))

