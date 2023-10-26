### xgboost on SACNZ - compare ngrams###
# start with a few genes

# Load libraries and functions ##########################################################################################################
source("methods/libs_fns.R")
library(tidymodels)
library(finetune)
#library(xgboost)i8
#library(vip)

# define functions #########################################################################################################
prepare_training_pco_m <- function(data, vars, class, d, m.list, mp, residualised=NULL) {
  # pull out our var_cols and class
  var_cols <- dplyr::select(data, all_of(vars))
  classes   <- data |> pull(class)
  # iterate over the var columns and distance matrices, and convert
  prepped <- pmap(list(var_cols, d, m.list), factor_to_pco_score, mp) |> compact() # removes empties
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
  var_cols <- dplyr::select(data,all_of(vars))
  classes   <- data |> pull(class)
  # iterate over the var columns and distance matrices, and convert
  prepped <- pmap(list(var_cols, d, m.list), factor_to_CAP_score, class = classes, k, mp, c) |> compact() # removes empties
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
load("../CAP_Data/data/list_of_6grams.RData")
load("../CAP_Data/data/list_of_8grams.RData")
load("../CAP_Data/data/list_of_10grams.RData")
load("../CAP_Data/data/list_of_12grams.RData")

# prepare data - pull out gene names and select animal only isolates
cgMLST_dat <- cgMLST |> 
  mutate(across(everything(), factor)) |> 
  mutate(Source = factor(Source, levels = c("Beef", "Poultry", "Sheep", "Human"))) |>
  mutate(Source = fct_recode(as.factor(Source), Cattle="Beef", Chicken="Poultry")) |> 
  filter(Source != "Human") |> droplevels()

# subset for testing
# dat <- cgMLST_dat |> select(c(LabID, Source, sample(2:1344, size=5)))
# d.6gram <- list_of_6grams[colnames(dat |> select(starts_with("CAMP")))]
# d.8gram <- list_of_8grams[colnames(dat |> select(starts_with("CAMP")))]
# d.10gram <- list_of_10grams[colnames(dat |> select(starts_with("CAMP")))]
# d.12gram <- list_of_12grams[colnames(dat |> select(starts_with("CAMP")))]

# all data
dat <- cgMLST_dat
d.6gram <- list_of_6grams
d.8gram <- list_of_8grams
d.10gram <- list_of_10grams
d.12gram <- list_of_12grams

# split into folds #########################################################################################################
set.seed(347)
folds <- vfold_cv(dat, strata = Source)

# score each fold #########################################################################################################
m.null <- list(NULL) |> rep(length(d.6gram))
names(m.null) <- names(d.6gram)

pco_6gram_folds <- folds
m <- list()
for(i in 1:10){
  x <- pco_6gram_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d.6gram, m.null, mp=95)$m
}
m.95 <- map(m, unlist) |> bind_rows() |> map(max)
m <- list()
for(i in 1:10){
  x <- pco_6gram_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d.6gram, m.null, mp=99)$m
}
m.99 <- map(m, unlist) |> bind_rows() |> map(max)
for(i in 1:10){
  x <- pco_6gram_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="pco", d=d.6gram, m.list=m.95, mp=NULL)
  pco_6gram_folds$splits[[i]]$data <- prepped.dat$dat
}
if((pco_6gram_folds$splits |> map(function(x) x$data |> ncol()) )|> unlist() |> table() |> as.data.frame() |> ncol() != 1){
  f <- pco_6gram_folds$splits |> map(function(x) x$data |> ncol()) |> which.min()
  cols <- pco_6gram_folds$splits[[f]]$data |> colnames()
  for(i in 1:10){
    pco_6gram_folds$splits[[i]]$data <- pco_6gram_folds$splits[[i]]$data |> select(all_of(cols))
  }
}


cap_6gram_folds <- folds # use m.99 as m.95 has a lot of single dimension variables
for(i in 1:10){
  x <- cap_6gram_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="cap", m.list=m.99, mp=NULL, d=d.6gram, k=2)
  cap_6gram_folds$splits[[i]]$data <- prepped.dat$dat
}
if((cap_6gram_folds$splits |> map(function(x) x$data |> ncol()) )|> unlist() |> table() |> as.data.frame() |> ncol() != 1){
  f <- cap_6gram_folds$splits |> map(function(x) x$data |> ncol()) |> which.min()
  cols <- cap_6gram_folds$splits[[f]]$data |> colnames()
  for(i in 1:10){
    cap_6gram_folds$splits[[i]]$data <- cap_6gram_folds$splits[[i]]$data |> select(all_of(cols))
  }
}

pco_8gram_folds <- folds
m <- list()
for(i in 1:10){
  x <- pco_8gram_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d.8gram, m.null, mp=95)$m
}
m.95 <- map(m, unlist) |> bind_rows() |> map(max)
m <- list()
for(i in 1:10){
  x <- pco_8gram_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d.8gram, m.null, mp=99)$m
}
m.99 <- map(m, unlist) |> bind_rows() |> map(max)
for(i in 1:10){
  x <- pco_8gram_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="pco", d=d.8gram, m.list=m.95, mp=NULL)
  pco_8gram_folds$splits[[i]]$data <- prepped.dat$dat
}
if((pco_8gram_folds$splits |> map(function(x) x$data |> ncol()) )|> unlist() |> table() |> as.data.frame() |> ncol() != 1){
  f <- pco_8gram_folds$splits |> map(function(x) x$data |> ncol()) |> which.min()
  cols <- pco_8gram_folds$splits[[f]]$data |> colnames()
  for(i in 1:10){
    pco_8gram_folds$splits[[i]]$data <- pco_8gram_folds$splits[[i]]$data |> select(all_of(cols))
  }
}


cap_8gram_folds <- folds # use m.99 as m.95 has a lot of single dimension variables
for(i in 1:10){
  x <- cap_8gram_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="cap", m.list=m.99, mp=NULL, d=d.8gram, k=2)
  cap_8gram_folds$splits[[i]]$data <- prepped.dat$dat
}
if((cap_8gram_folds$splits |> map(function(x) x$data |> ncol()) )|> unlist() |> table() |> as.data.frame() |> ncol() != 1){
  f <- cap_8gram_folds$splits |> map(function(x) x$data |> ncol()) |> which.min()
  cols <- cap_8gram_folds$splits[[f]]$data |> colnames()
  for(i in 1:10){
    cap_8gram_folds$splits[[i]]$data <- cap_8gram_folds$splits[[i]]$data |> select(all_of(cols))
  }
}

pco_10gram_folds <- folds
m <- list()
for(i in 1:10){
  x <- pco_10gram_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d.10gram, m.null, mp=95)$m
}
m.95 <- map(m, unlist) |> bind_rows() |> map(max)
m <- list()
for(i in 1:10){
  x <- pco_10gram_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d.10gram, m.null, mp=99)$m
}
m.99 <- map(m, unlist) |> bind_rows() |> map(max)
for(i in 1:10){
  x <- pco_10gram_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="pco", d=d.10gram, m.list=m.95, mp=NULL)
  pco_10gram_folds$splits[[i]]$data <- prepped.dat$dat
}
if((pco_10gram_folds$splits |> map(function(x) x$data |> ncol()) )|> unlist() |> table() |> as.data.frame() |> ncol() != 1){
  f <- pco_10gram_folds$splits |> map(function(x) x$data |> ncol()) |> which.min()
  cols <- pco_10gram_folds$splits[[f]]$data |> colnames()
  for(i in 1:10){
    pco_10gram_folds$splits[[i]]$data <- pco_10gram_folds$splits[[i]]$data |> select(all_of(cols))
  }
}


cap_10gram_folds <- folds # use m.99 as m.95 has a lot of single dimension variables
for(i in 1:10){
  x <- cap_10gram_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="cap", m.list=m.99, mp=NULL, d=d.10gram, k=2)
  cap_10gram_folds$splits[[i]]$data <- prepped.dat$dat
}
if((cap_10gram_folds$splits |> map(function(x) x$data |> ncol()) )|> unlist() |> table() |> as.data.frame() |> ncol() != 1){
  f <- cap_10gram_folds$splits |> map(function(x) x$data |> ncol()) |> which.min()
  cols <- cap_10gram_folds$splits[[f]]$data |> colnames()
  for(i in 1:10){
    cap_10gram_folds$splits[[i]]$data <- cap_10gram_folds$splits[[i]]$data |> select(all_of(cols))
  }
}

pco_12gram_folds <- folds
m <- list()
for(i in 1:10){
  x <- pco_12gram_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d.12gram, m.null, mp=95)$m
}
m.95 <- map(m, unlist) |> bind_rows() |> map(max)
m <- list()
for(i in 1:10){
  x <- pco_12gram_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d.12gram, m.null, mp=99)$m
}
m.99 <- map(m, unlist) |> bind_rows() |> map(max)
for(i in 1:10){
  x <- pco_12gram_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="pco", d=d.12gram, m.list=m.95, mp=NULL)
  pco_12gram_folds$splits[[i]]$data <- prepped.dat$dat
}
if((pco_12gram_folds$splits |> map(function(x) x$data |> ncol()) )|> unlist() |> table() |> as.data.frame() |> ncol() != 1){
  f <- pco_12gram_folds$splits |> map(function(x) x$data |> ncol()) |> which.min()
  cols <- pco_12gram_folds$splits[[f]]$data |> colnames()
  for(i in 1:10){
    pco_12gram_folds$splits[[i]]$data <- pco_12gram_folds$splits[[i]]$data |> select(all_of(cols))
  }
}


cap_12gram_folds <- folds # use m.99 as m.95 has a lot of single dimension variables
for(i in 1:10){
  x <- cap_12gram_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="cap", m.list=m.99, mp=NULL, d=d.12gram, k=2)
  cap_12gram_folds$splits[[i]]$data <- prepped.dat$dat
}
if((cap_12gram_folds$splits |> map(function(x) x$data |> ncol()) )|> unlist() |> table() |> as.data.frame() |> ncol() != 1){
  f <- cap_12gram_folds$splits |> map(function(x) x$data |> ncol()) |> which.min()
  cols <- cap_12gram_folds$splits[[f]]$data |> colnames()
  for(i in 1:10){
    cap_12gram_folds$splits[[i]]$data <- cap_12gram_folds$splits[[i]]$data |> select(all_of(cols))
  }
}




# define xgboost model, recipe and workflow #########################################################################################################
xg_spec <- 
  boost_tree(trees = tune(), #default 15
             min_n = tune(), #default 1
             mtry = tune(), #default 
             #tree_depth = tune(), #default 6           
             #learn_rate = 0.01
             learn_rate = 0.1) |>
  set_engine("xgboost") |> 
  set_mode("classification")

xg_recipe.pco.6g <- recipe(Source ~ ., data = pco_6gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.pco.6g <- workflow(xg_recipe.pco.6g, xg_spec)

xg_recipe.cap.6g <- recipe(Source ~ ., data = cap_6gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.cap.6g <- workflow(xg_recipe.cap.6g, xg_spec)

xg_recipe.pco.8g <- recipe(Source ~ ., data = pco_8gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.pco.8g <- workflow(xg_recipe.pco.8g, xg_spec)

xg_recipe.cap.8g <- recipe(Source ~ ., data = cap_8gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.cap.8g <- workflow(xg_recipe.cap.8g, xg_spec)

xg_recipe.pco.10g <- recipe(Source ~ ., data = pco_10gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.pco.10g <- workflow(xg_recipe.pco.10g, xg_spec)

xg_recipe.cap.10g <- recipe(Source ~ ., data = cap_10gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.cap.10g <- workflow(xg_recipe.cap.10g, xg_spec)

xg_recipe.pco.12g <- recipe(Source ~ ., data = pco_12gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.pco.12g <- workflow(xg_recipe.pco.12g, xg_spec)

xg_recipe.cap.12g <- recipe(Source ~ ., data = cap_12gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.cap.12g <- workflow(xg_recipe.cap.12g, xg_spec)


# tune xgboost parameters ##############################################################################################################
doParallel::registerDoParallel()
set.seed(319)

xg_resamples.pco.6g <- tune_race_anova(
  xg_workflow.pco.6g,
  resamples = pco_6gram_folds,
  grid = 15,
  control = control_race(verbose_elim = TRUE)
)

xg_resamples.cap.6g <- tune_race_anova(
  xg_workflow.cap.6g,
  resamples = cap_6gram_folds,
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

xg_resamples.pco.10g <- tune_race_anova(
  xg_workflow.pco.10g,
  resamples = pco_10gram_folds,
  grid = 15,
  control = control_race(verbose_elim = TRUE)
)

xg_resamples.cap.10g <- tune_race_anova(
  xg_workflow.cap.10g,
  resamples = cap_10gram_folds,
  grid = 15,
  control = control_race(verbose_elim = TRUE)
)

xg_resamples.pco.12g <- tune_race_anova(
  xg_workflow.pco.12g,
  resamples = pco_12gram_folds,
  grid = 15,
  control = control_race(verbose_elim = TRUE)
)

xg_resamples.cap.12g <- tune_race_anova(
  xg_workflow.cap.12g,
  resamples = cap_12gram_folds,
  grid = 15,
  control = control_race(verbose_elim = TRUE)
)


# find best parameters ##############################################################################################################
best_pco.6g <- select_best(xg_resamples.pco.6g, "accuracy")
best_cap.6g <- select_best(xg_resamples.cap.6g, "accuracy")
best_pco.8g <- select_best(xg_resamples.pco.8g, "accuracy")
best_cap.8g <- select_best(xg_resamples.cap.8g, "accuracy")
best_pco.10g <- select_best(xg_resamples.pco.10g, "accuracy")
best_cap.10g <- select_best(xg_resamples.cap.10g, "accuracy")
best_pco.12g <- select_best(xg_resamples.pco.12g, "accuracy")
best_cap.12g <- select_best(xg_resamples.cap.12g, "accuracy")


# split into folds ##############################################################################################################
set.seed(795)
folds <- vfold_cv(dat, strata = Source)

# score each fold ##############################################################################################################
m.null <- list(NULL) |> rep(length(d.6gram))
names(m.null) <- names(d.6gram)

pco_6gram_folds <- folds
m <- list()
for(i in 1:10){
  x <- pco_6gram_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d.6gram, m.null, mp=95)$m
}
m.95 <- map(m, unlist) |> bind_rows() |> map(max)
m <- list()
for(i in 1:10){
  x <- pco_6gram_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d.6gram, m.null, mp=99)$m
}
m.99 <- map(m, unlist) |> bind_rows() |> map(max)
for(i in 1:10){
  x <- pco_6gram_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="pco", d=d.6gram, m.list=m.95, mp=NULL)
  pco_6gram_folds$splits[[i]]$data <- prepped.dat$dat
}
if((pco_6gram_folds$splits |> map(function(x) x$data |> ncol()) )|> unlist() |> table() |> as.data.frame() |> ncol() != 1){
  f <- pco_6gram_folds$splits |> map(function(x) x$data |> ncol()) |> which.min()
  cols <- pco_6gram_folds$splits[[f]]$data |> colnames()
  for(i in 1:10){
    pco_6gram_folds$splits[[i]]$data <- pco_6gram_folds$splits[[i]]$data |> select(all_of(cols))
  }
}

cap_6gram_folds <- folds # use m.99 as m.95 has a lot of single dimension variables
for(i in 1:10){
  x <- cap_6gram_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="cap", m.list=m.99, mp=NULL, d=d.6gram)
  cap_6gram_folds$splits[[i]]$data <- prepped.dat$dat
}
cap_6gram_folds$splits |> map(function(x) x$data |> ncol()) 
cap_6gram_folds$splits |> map(function(x) x$data |> ncol()) 
if((cap_6gram_folds$splits |> map(function(x) x$data |> ncol()) )|> unlist() |> table() |> as.data.frame() |> ncol() != 1){
  f <- cap_6gram_folds$splits |> map(function(x) x$data |> ncol()) |> which.min()
  cols <- cap_6gram_folds$splits[[f]]$data |> colnames()
  for(i in 1:10){
    cap_6gram_folds$splits[[i]]$data <- cap_6gram_folds$splits[[i]]$data |> select(all_of(cols))
  }
}

pco_8gram_folds <- folds
m <- list()
for(i in 1:10){
  x <- pco_8gram_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d.8gram, m.null, mp=95)$m
}
m.95 <- map(m, unlist) |> bind_rows() |> map(max)
m <- list()
for(i in 1:10){
  x <- pco_8gram_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d.8gram, m.null, mp=99)$m
}
m.99 <- map(m, unlist) |> bind_rows() |> map(max)
for(i in 1:10){
  x <- pco_8gram_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="pco", d=d.8gram, m.list=m.95, mp=NULL)
  pco_8gram_folds$splits[[i]]$data <- prepped.dat$dat
}
if((pco_8gram_folds$splits |> map(function(x) x$data |> ncol()) )|> unlist() |> table() |> as.data.frame() |> ncol() != 1){
  f <- pco_8gram_folds$splits |> map(function(x) x$data |> ncol()) |> which.min()
  cols <- pco_8gram_folds$splits[[f]]$data |> colnames()
  for(i in 1:10){
    pco_8gram_folds$splits[[i]]$data <- pco_8gram_folds$splits[[i]]$data |> select(all_of(cols))
  }
}

cap_8gram_folds <- folds # use m.99 as m.95 has a lot of single dimension variables
for(i in 1:10){
  x <- cap_8gram_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="cap", m.list=m.99, mp=NULL, d=d.8gram)
  cap_8gram_folds$splits[[i]]$data <- prepped.dat$dat
}
cap_8gram_folds$splits |> map(function(x) x$data |> ncol()) 
cap_8gram_folds$splits |> map(function(x) x$data |> ncol()) 
if((cap_8gram_folds$splits |> map(function(x) x$data |> ncol()) )|> unlist() |> table() |> as.data.frame() |> ncol() != 1){
  f <- cap_8gram_folds$splits |> map(function(x) x$data |> ncol()) |> which.min()
  cols <- cap_8gram_folds$splits[[f]]$data |> colnames()
  for(i in 1:10){
    cap_8gram_folds$splits[[i]]$data <- cap_8gram_folds$splits[[i]]$data |> select(all_of(cols))
  }
}

pco_10gram_folds <- folds
m <- list()
for(i in 1:10){
  x <- pco_10gram_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d.10gram, m.null, mp=95)$m
}
m.95 <- map(m, unlist) |> bind_rows() |> map(max)
m <- list()
for(i in 1:10){
  x <- pco_10gram_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d.10gram, m.null, mp=99)$m
}
m.99 <- map(m, unlist) |> bind_rows() |> map(max)
for(i in 1:10){
  x <- pco_10gram_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="pco", d=d.10gram, m.list=m.95, mp=NULL)
  pco_10gram_folds$splits[[i]]$data <- prepped.dat$dat
}
if((pco_10gram_folds$splits |> map(function(x) x$data |> ncol()) )|> unlist() |> table() |> as.data.frame() |> ncol() != 1){
  f <- pco_10gram_folds$splits |> map(function(x) x$data |> ncol()) |> which.min()
  cols <- pco_10gram_folds$splits[[f]]$data |> colnames()
  for(i in 1:10){
    pco_10gram_folds$splits[[i]]$data <- pco_10gram_folds$splits[[i]]$data |> select(all_of(cols))
  }
}

cap_10gram_folds <- folds # use m.99 as m.95 has a lot of single dimension variables
for(i in 1:10){
  x <- cap_10gram_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="cap", m.list=m.99, mp=NULL, d=d.10gram)
  cap_10gram_folds$splits[[i]]$data <- prepped.dat$dat
}
cap_10gram_folds$splits |> map(function(x) x$data |> ncol()) 
cap_10gram_folds$splits |> map(function(x) x$data |> ncol()) 
if((cap_10gram_folds$splits |> map(function(x) x$data |> ncol()) )|> unlist() |> table() |> as.data.frame() |> ncol() != 1){
  f <- cap_10gram_folds$splits |> map(function(x) x$data |> ncol()) |> which.min()
  cols <- cap_10gram_folds$splits[[f]]$data |> colnames()
  for(i in 1:10){
    cap_10gram_folds$splits[[i]]$data <- cap_10gram_folds$splits[[i]]$data |> select(all_of(cols))
  }
}


pco_12gram_folds <- folds
m <- list()
for(i in 1:10){
  x <- pco_12gram_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d.12gram, m.null, mp=95)$m
}
m.95 <- map(m, unlist) |> bind_rows() |> map(max)
m <- list()
for(i in 1:10){
  x <- pco_12gram_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d.12gram, m.null, mp=99)$m
}
m.99 <- map(m, unlist) |> bind_rows() |> map(max)
for(i in 1:10){
  x <- pco_12gram_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="pco", d=d.12gram, m.list=m.95, mp=NULL)
  pco_12gram_folds$splits[[i]]$data <- prepped.dat$dat
}
if((pco_12gram_folds$splits |> map(function(x) x$data |> ncol()) )|> unlist() |> table() |> as.data.frame() |> ncol() != 1){
  f <- pco_12gram_folds$splits |> map(function(x) x$data |> ncol()) |> which.min()
  cols <- pco_12gram_folds$splits[[f]]$data |> colnames()
  for(i in 1:10){
    pco_12gram_folds$splits[[i]]$data <- pco_12gram_folds$splits[[i]]$data |> select(all_of(cols))
  }
}

cap_12gram_folds <- folds # use m.99 as m.95 has a lot of single dimension variables
for(i in 1:10){
  x <- cap_12gram_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="cap", m.list=m.99, mp=NULL, d=d.12gram)
  cap_12gram_folds$splits[[i]]$data <- prepped.dat$dat
}
cap_12gram_folds$splits |> map(function(x) x$data |> ncol()) 
cap_12gram_folds$splits |> map(function(x) x$data |> ncol()) 
if((cap_12gram_folds$splits |> map(function(x) x$data |> ncol()) )|> unlist() |> table() |> as.data.frame() |> ncol() != 1){
  f <- cap_12gram_folds$splits |> map(function(x) x$data |> ncol()) |> which.min()
  cols <- cap_12gram_folds$splits[[f]]$data |> colnames()
  for(i in 1:10){
    cap_12gram_folds$splits[[i]]$data <- cap_12gram_folds$splits[[i]]$data |> select(all_of(cols))
  }
}


# update model specs ############################################################################################################
# need to update the recipe as well as the columns might be different in the new data folds
xg_recipe.pco.6g <- recipe(Source ~ ., data = pco_6gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.pco.6g <- workflow(xg_recipe.pco.6g, xg_spec)
xg_final_workflow.pco.6g <- finalize_workflow(xg_workflow.pco.6g, best_pco.6g)

xg_recipe.cap.6g <- recipe(Source ~ ., data = cap_6gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.cap.6g <- workflow(xg_recipe.cap.6g, xg_spec)
xg_final_workflow.cap.6g <- finalize_workflow(xg_workflow.cap.6g, best_cap.6g)

xg_recipe.pco.8g <- recipe(Source ~ ., data = pco_8gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.pco.8g <- workflow(xg_recipe.pco.8g, xg_spec)
xg_final_workflow.pco.8g <- finalize_workflow(xg_workflow.pco.8g, best_pco.8g)

xg_recipe.cap.8g <- recipe(Source ~ ., data = cap_8gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.cap.8g <- workflow(xg_recipe.cap.8g, xg_spec)
xg_final_workflow.cap.8g <- finalize_workflow(xg_workflow.cap.8g, best_cap.8g)

xg_recipe.pco.10g <- recipe(Source ~ ., data = pco_10gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.pco.10g <- workflow(xg_recipe.pco.10g, xg_spec)
xg_final_workflow.pco.10g <- finalize_workflow(xg_workflow.pco.10g, best_pco.10g)

xg_recipe.cap.10g <- recipe(Source ~ ., data = cap_10gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.cap.10g <- workflow(xg_recipe.cap.10g, xg_spec)
xg_final_workflow.cap.10g <- finalize_workflow(xg_workflow.cap.10g, best_cap.10g)

xg_recipe.pco.12g <- recipe(Source ~ ., data = pco_12gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.pco.12g <- workflow(xg_recipe.pco.12g, xg_spec)
xg_final_workflow.pco.12g <- finalize_workflow(xg_workflow.pco.12g, best_pco.12g)

xg_recipe.cap.12g <- recipe(Source ~ ., data = cap_12gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
xg_workflow.cap.12g <- workflow(xg_recipe.cap.12g, xg_spec)
xg_final_workflow.cap.12g <- finalize_workflow(xg_workflow.cap.12g, best_cap.12g)


# run xgboost resamples ##############################################################################################################
set.seed(273)
xg_rs.pco.6g <- fit_resamples(
  xg_final_workflow.pco.6g,
  resamples = pco_6gram_folds,
  control = control_resamples(save_pred = TRUE)
)

xg_rs.cap.6g <- fit_resamples(
  xg_final_workflow.cap.6g,
  resamples = cap_6gram_folds,
  control = control_resamples(save_pred = TRUE)
)

xg_rs.pco.8g <- fit_resamples(
  xg_final_workflow.pco.8g,
  resamples = pco_8gram_folds,
  control = control_resamples(save_pred = TRUE)
)

xg_rs.cap.8g <- fit_resamples(
  xg_final_workflow.cap.8g,
  resamples = cap_8gram_folds,
  control = control_resamples(save_pred = TRUE)
)

xg_rs.pco.10g <- fit_resamples(
  xg_final_workflow.pco.10g,
  resamples = pco_10gram_folds,
  control = control_resamples(save_pred = TRUE)
)

xg_rs.cap.10g <- fit_resamples(
  xg_final_workflow.cap.10g,
  resamples = cap_10gram_folds,
  control = control_resamples(save_pred = TRUE)
)

xg_rs.pco.12g <- fit_resamples(
  xg_final_workflow.pco.12g,
  resamples = pco_12gram_folds,
  control = control_resamples(save_pred = TRUE)
)

xg_rs.cap.12g <- fit_resamples(
  xg_final_workflow.cap.12g,
  resamples = cap_12gram_folds,
  control = control_resamples(save_pred = TRUE)
)

# define rf model, recipe and workflow ##############################################################################################################
rf_spec <- 
  rand_forest(trees = 500)  |>  
  set_engine("ranger", respect.unordered.factors = TRUE) |> 
  set_mode("classification")

rf_recipe.pco.6g <- recipe(Source ~ ., data = pco_6gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
rf_workflow.pco.6g <- workflow(rf_recipe.pco.6g, rf_spec)

rf_recipe.cap.6g <- recipe(Source ~ ., data = cap_6gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
rf_workflow.cap.6g <- workflow(rf_recipe.cap.6g, rf_spec)

rf_recipe.pco.8g <- recipe(Source ~ ., data = pco_8gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
rf_workflow.pco.8g <- workflow(rf_recipe.pco.8g, rf_spec)

rf_recipe.cap.8g <- recipe(Source ~ ., data = cap_8gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
rf_workflow.cap.8g <- workflow(rf_recipe.cap.8g, rf_spec)

rf_recipe.pco.10g <- recipe(Source ~ ., data = pco_10gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
rf_workflow.pco.10g <- workflow(rf_recipe.pco.10g, rf_spec)

rf_recipe.cap.10g <- recipe(Source ~ ., data = cap_10gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
rf_workflow.cap.10g <- workflow(rf_recipe.cap.10g, rf_spec)

rf_recipe.pco.12g <- recipe(Source ~ ., data = pco_12gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
rf_workflow.pco.12g <- workflow(rf_recipe.pco.12g, rf_spec)

rf_recipe.cap.12g <- recipe(Source ~ ., data = cap_12gram_folds$splits[[1]]$data) |> 
  update_role(LabID, new_role = "id")
rf_workflow.cap.12g <- workflow(rf_recipe.cap.12g, rf_spec)


# run rf resamples ##############################################################################################################
set.seed(288)
rf_rs.pco.6g <- fit_resamples(
  rf_workflow.pco.6g,
  resamples = pco_6gram_folds,
  control = control_resamples(save_pred = TRUE)
)

rf_rs.cap.6g <- fit_resamples(
  rf_workflow.cap.6g,
  resamples = cap_6gram_folds,
  control = control_resamples(save_pred = TRUE)
)

rf_rs.pco.8g <- fit_resamples(
  rf_workflow.pco.8g,
  resamples = pco_8gram_folds,
  control = control_resamples(save_pred = TRUE)
)

rf_rs.cap.8g <- fit_resamples(
  rf_workflow.cap.8g,
  resamples = cap_8gram_folds,
  control = control_resamples(save_pred = TRUE)
)

rf_rs.pco.10g <- fit_resamples(
  rf_workflow.pco.10g,
  resamples = pco_10gram_folds,
  control = control_resamples(save_pred = TRUE)
)

rf_rs.cap.10g <- fit_resamples(
  rf_workflow.cap.10g,
  resamples = cap_10gram_folds,
  control = control_resamples(save_pred = TRUE)
)

rf_rs.pco.12g <- fit_resamples(
  rf_workflow.pco.12g,
  resamples = pco_12gram_folds,
  control = control_resamples(save_pred = TRUE)
)

rf_rs.cap.12g <- fit_resamples(
  rf_workflow.cap.12g,
  resamples = cap_12gram_folds,
  control = control_resamples(save_pred = TRUE)
)

ngram_metrics <- list(pco_6g=list(rf=rf_rs.pco.6g, xg=xg_rs.pco.6g),
                      cap_6g=list(rf=rf_rs.cap.6g,xg=xg_rs.cap.6g),
                      pco_8g=list(rf=rf_rs.pco.8g, xg=xg_rs.pco.8g),
                      cap_8g=list(rf=rf_rs.cap.8g,xg=xg_rs.cap.8g),
                      pco_10g=list(rf=rf_rs.pco.10g, xg=xg_rs.pco.10g),
                      cap_10g=list(rf=rf_rs.cap.10g,xg=xg_rs.cap.10g),
                      pco_12g=list(rf=rf_rs.pco.12g, xg=xg_rs.pco.12g),
                      cap_12g=list(rf=rf_rs.cap.12g,xg=xg_rs.cap.12g))
#save(ngram_metrics, file = "ngram_metrics_SACNZ_v2.Rdata")
#load("ngram_metrics_SACNZ_v2.Rdata")
#save(ngram_metrics, file = "ngram_metrics_SACNZ.Rdata")
#load("ngram_metrics_SACNZ.Rdata")

# plot ##############################################################################################################

p <- ngram_metrics |> map_df(function(x) x |> map_df(collect_metrics, .id = "model"), .id = "method" ) |>
  filter(.metric == "accuracy") |>
  mutate(dist = factor(case_when(method |> str_detect("6g") ~ "6g",
                                 method |> str_detect("8g") ~ "8g",
                                 method |> str_detect("10g") ~ "10g",
                                 method |> str_detect("12g") ~ "12g"), ordered=TRUE, levels=c("6g","8g","10g","12g")),
         method2 = factor(method),
         method = factor(str_sub(method2, 1,3), ordered=TRUE, levels=c("pco", "cap"))) |>
  select(method2, method, model, mean, std_err, dist)

p |>
  ggplot(aes(x=method, y=mean)) +
  geom_errorbar(aes(col=dist, group=as.factor(dist), ymin=mean-std_err, ymax=mean+std_err), width=0.2, position = position_dodge(0.5)) +
  geom_point(aes(col=dist, group=as.factor(dist)), size=3, stroke=1.5, position = position_dodge(0.5)) +
  #ylim(0,1) +
  #  scale_colour_manual(values = c("#115896", "#4C4C53", "#BA2F00")) +
  guides(col = guide_legend("Distance measure")) +
  labs(x = "Method of encoding", y = "Classification success", title = "xgboost and rf results, SACNZ sample data, 10 cv folds") +
  theme_bw() +
  facet_wrap(~model, labeller = labeller(model=c("rf" = "Random Forest", "xg" = "xgboost")))

p2 <- bind_rows(
  rf = ngram_metrics |> map_df(function(x) x$rf |> pull(.metrics) |> map_df(function(x) x |> filter(.metric=="accuracy") |> select(.estimate), .id="fold"), .id = "method" ),
  xg = ngram_metrics |> map_df(function(x) x$xg |> pull(.metrics) |> map_df(function(x) x |> filter(.metric=="accuracy") |> select(.estimate), .id="fold"), .id = "method" ), .id="model") |>
  mutate(dist = factor(case_when(method |> str_detect("6g") ~ "6g",
                                 method |> str_detect("8g") ~ "8g",
                                 method |> str_detect("10g") ~ "10g",
                                 method |> str_detect("12g") ~ "12g"), ordered=TRUE, levels=c("6g","8g","10g","12g")),
         method2 = factor(method),
         method = factor(str_sub(method2, 1,3), ordered=TRUE, levels=c("pco", "cap")))

p2 |>
  ggplot(aes(x=method, y=.estimate)) +
  geom_point(aes(col=dist, group=as.factor(dist)), size=3, stroke=1.5, position = position_dodge(0.5)) +
  stat_summary(aes(group=as.factor(dist), col=dist), geom = "point", fun = "mean", size = 5, shape = 4, position = position_dodge(0.5)) +
  ylim(0,1) +
  #  scale_colour_manual(values = c("#115896", "#4C4C53", "#BA2F00")) +
  guides(col = guide_legend("Distance measure"),
         shape = guide_legend("Distance measure")) +
  labs(x = "Method of encoding", y = "Classification success", title = "xgboost and rf results, SACNZ sample data, 10 cv folds") +
  theme_bw() +
  facet_wrap(~model, labeller = labeller(model=c("rf" = "Random Forest", "xg" = "xgboost")))

p3 <- bind_rows(
        rf = ngram_metrics |> map_df(function(x) x$rf$.predictions |> reduce(rbind) |> select(.pred_class, Source) |> group_by(Source) |>
                                 summarise(N=n(), n=sum(.pred_class==Source), p=n/N), .id="method"),
        xg = ngram_metrics |> map_df(function(x) x$xg$.predictions |> reduce(rbind) |> select(.pred_class, Source) |> group_by(Source) |>
                                 summarise(N=n(), n=sum(.pred_class==Source), p=n/N), .id="method"), .id = "model") |>
      mutate(dist = factor(case_when(method |> str_detect("6g") ~ "6g",
                                 method |> str_detect("8g") ~ "8g",
                                 method |> str_detect("10g") ~ "10g",
                                 method |> str_detect("12g") ~ "12g"), ordered=TRUE, levels=c("6g","8g","10g","12g")),
             method2 = factor(method),
             method = factor(str_sub(method2, 1,3), ordered=TRUE, levels=c("pco", "cap")))

p3 |>
  ggplot(aes(x=method, y=p)) +
  geom_point(aes(col=dist, shape=Source, group=dist), size=3, stroke=1.5, position = position_dodge(0.5)) +
  ylim(0,1) +
  #scale_colour_manual(values = c("#115896", "#4C4C53", "#BA2F00")) +
  scale_shape_manual(values = c("Sheep"=1, "Cattle"=2, "Chicken"=0)) +
  guides(col = guide_legend("ngram"), shape = guide_legend("Source")) +
  labs(x = "Method of encoding", y = "Classification success", title = "xgboost and rf results, SACNZ sample data, 10 cv folds") +
  theme_bw() +
  facet_wrap(~model, labeller = labeller(model=c("rf" = "Random Forest", "xg" = "xgboost")))

