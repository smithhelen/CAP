### xgboost on SACNZ ###
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

# prepare data - pull out gene names and select animal only isolates
cgMLST_dat <- cgMLST |> filter(Source != "Human") |> droplevels()

# distance information
load("../CAP_Data/data/list_of_distance_matrices.RData") #  distance matrices of hamming distances among unique alleles for each gene

# all data
dat <- cgMLST_dat

# split into folds #########################################################################################################
set.seed(2471)
folds <- vfold_cv(dat, strata = Source)

# score each fold #########################################################################################################
ca0_folds <- folds
for(i in 1:10){
  x <- ca0_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="ca0", d=NULL)
  ca0_folds$splits[[i]]$data <- prepped.dat$dat
}

pco_folds <- folds
# find the maximum number of pco axes for each gene

# define list of m=null so that the number of axes will be determined by propG/mp
m.null <- list(NULL) |> rep(length(d))
names(m.null) <- names(d)

# create empty list for mp=95%
m <- list()
for(i in 1:10){
  x <- pco_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d, m.null, mp=95)$m
}
m.95 <- map(m, unlist) |> bind_rows() |> map(max)

# create empty list for mp=99%
m <- list()
for(i in 1:10){
  x <- pco_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d, m.null, mp=99)$m
}
m.99 <- map(m, unlist) |> bind_rows() |> map(max)

# now use m (number of axes) rather than mp (prop var)
for(i in 1:10){
  x <- pco_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="pco", d=d, m.list=m.95, mp=NULL)
  pco_folds$splits[[i]]$data <- prepped.dat$dat
}
#pco_folds$splits |> map(function(x) x$data |> ncol()) 
#pco_folds$splits |> map(function(x) x$data |> colnames()) 


cap_folds <- folds # use m.99 as m.95 has a lot of single dimension variables
for(i in 1:10){
  x <- cap_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="cap", m.list=m.99, mp=NULL, d=d)
  cap_folds$splits[[i]]$data <- prepped.dat$dat
}
# cap_folds <- folds # m.95
# for(i in 1:10){
#   x <- cap_folds$splits[[i]]
#   prepped.dat <- score_fn_m(x, method="cap", m.list=m.95, mp=NULL, d=d)
#   cap_folds$splits[[i]]$data <- prepped.dat$dat
# }
#cap_folds$splits |> map(function(x) x$data |> ncol()) 
#cap_folds$splits |> map(function(x) x$data |> colnames()) 



# define xgboost model, recipe and workflow #########################################################################################################
xg_spec <- 
  boost_tree(trees = tune(), #default 15
             min_n = tune(), #default 1
             mtry = tune(), #default 
             #tree_depth = tune(), #default 6           
             learn_rate = 0.01) |>
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


# tune xgboost parameters ##############################################################################################################
doParallel::registerDoParallel()
set.seed(457)
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



# find best parameters ##############################################################################################################
#collect_metrics(xg_resamples.ca0)
#plot_race(xg_resamples.ca0)
#show_best(xg_resamples.ca0, "accuracy")
best_ca0 <- select_best(xg_resamples.ca0, "accuracy")

#collect_metrics(xg_resamples.pco)
#plot_race(xg_resamples.pco)
#show_best(xg_resamples.pco, "accuracy")
best_pco <- select_best(xg_resamples.pco, "accuracy")

#collect_metrics(xg_resamples.cap)
#plot_race(xg_resamples.cap)
#show_best(xg_resamples.cap, "accuracy")
best_cap <- select_best(xg_resamples.cap, "accuracy")



# split into folds ##############################################################################################################
set.seed(377)
folds <- vfold_cv(dat, strata = Source)

# score each fold ##############################################################################################################
ca0_folds <- folds
for(i in 1:10){
  x <- ca0_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="ca0", d=NULL)
  ca0_folds$splits[[i]]$data <- prepped.dat$dat
}

pco_folds <- folds
# find the maximum number of pco axes for each gene

# define list of m=null so that the number of axes will be determined by propG/mp
m.null <- list(NULL) |> rep(length(d))
names(m.null) <- names(d)

# create empty list for mp=95%
m <- list()
for(i in 1:10){
  x <- pco_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d, m.null, mp=95)$m
}
m.95 <- map(m, unlist) |> bind_rows() |> map(max)

# create empty list for mp=99%
m <- list()
for(i in 1:10){
  x <- pco_folds$splits[[i]]
  m[[i]] <- score_fn_m(x, method="pco", d=d, m.null, mp=99)$m
}
m.99 <- map(m, unlist) |> bind_rows() |> map(max)

# now use m (number of axes) rather than mp (prop var)
for(i in 1:10){
  x <- pco_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="pco", d=d, m.list=m.95, mp=NULL)
  pco_folds$splits[[i]]$data <- prepped.dat$dat
}
#pco_folds$splits |> map(function(x) x$data |> ncol()) 
#pco_folds$splits |> map(function(x) x$data |> colnames()) 


cap_folds <- folds # use m.99 as m.95 has a lot of single dimension variables
for(i in 1:10){
  x <- cap_folds$splits[[i]]
  prepped.dat <- score_fn_m(x, method="cap", m.list=m.99, mp=NULL, d=d)
  cap_folds$splits[[i]]$data <- prepped.dat$dat
}
# cap_folds <- folds # m.95
# for(i in 1:10){
#   x <- cap_folds$splits[[i]]
#   prepped.dat <- score_fn_m(x, method="cap", m.list=m.95, mp=NULL, d=d)
#   cap_folds$splits[[i]]$data <- prepped.dat$dat
# }
#cap_folds$splits |> map(function(x) x$data |> ncol()) 
#cap_folds$splits |> map(function(x) x$data |> colnames()) 




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


xg_workflows <- list(ca0=xg_final_workflow.ca0,
                     pco=xg_final_workflow.pco,
                     cap=xg_final_workflow.cap)
# save(all_workflows, file = "final_workflows_SACNZ.Rdata")

# run xgboost resamples ##############################################################################################################
set.seed(767)
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

xgb_metrics <- list(ca0=xg_rs.ca0, 
                    pco=xg_rs.pco, 
                    cap=xg_rs.cap)
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


# run rf resamples ##############################################################################################################

set.seed(787)
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

rf_metrics <- list(ca0=rf_rs.ca0,
                   pco=rf_rs.pco,
                   cap=rf_rs.cap
                   )
#save(rf_metrics, file = "rf_metrics_SACNZ.Rdata")

all_metrics <- list(ca0=list(rf=rf_rs.ca0, xg=xg_rs.ca0),
                    pco=list(rf=rf_rs.pco, xg=xg_rs.pco),
                    cap=list(rf=rf_rs.cap, xg=xg_rs.cap)
                    )
#save(all_metrics, file = "all_metrics_SACNZ.Rdata")


# plot ##############################################################################################################

p <- all_metrics |> map_df(function(x) x |> map_df(collect_metrics, .id = "model"), .id = "method" ) |> 
  filter(.metric == "accuracy") |> 
  mutate(pch = case_when(model == "rf" ~ 21,
                         model == "xg" ~ 24),
         method2 = factor(method)) 

p |> 
  ggplot() +
  geom_point(aes(x=as.numeric(method), y=mean, col=method, shape=pch), size=3, stroke=1.5) +
  scale_shape_identity(guide = "legend", breaks=c(21,24), labels=c("random forest","xgboost")) +
  scale_colour_manual(values = c("#115896", "#4C4C53", "#BA2F00")) +
  guides(col = guide_legend("Encoding method"),
         shape = guide_legend("Model", override.aes = list(shape=c(1,2), color="black"))) +
  labs(x = NULL, y = "Classification success", title = "xgboost and rf results, midwest survey data, 10 cv folds") +
  theme_bw() + scale_x_discrete(limits = c("CA0", "CAP", "PCO"))

p |> 
  ggplot(aes(x=as.numeric(method), y=mean)) +
  geom_errorbar(aes(col=method, ymin=,ymax=), width=0.2, position = position_dodge(0.3)) +
  geom_point(aes(col=method, shape=pch), size=3, stroke=1.5, position = position_dodge(0.3)) +
  scale_shape_identity(guide = "legend", breaks=c(21,24), labels=c("random forest","xgboost")) +
  scale_colour_manual(values = c("#115896", "#4C4C53", "#BA2F00")) +
  guides(col = guide_legend("Encoding method"),
         shape = guide_legend("Model", override.aes = list(shape=c(1,2), color="black"))) +
  labs(x = NULL, y = "Classification success", title = "xgboost and rf results, midwest survey data, 10 cv folds") +
  theme_bw() + scale_x_discrete(limits = c("CA0", "CAP", "PCO"))

p |> 
  ggplot(aes(x=as.numeric(method), y=mean)) +
  geom_pointrange(aes(ymin=,ymax=,col=method, shape=pch), size=3, stroke=1.5, position = position_dodge(0.3)) +
  scale_shape_identity(guide = "legend", breaks=c(21,24), labels=c("random forest","xgboost")) +
  scale_colour_manual(values = c("#115896", "#4C4C53", "#BA2F00")) +
  guides(col = guide_legend("Encoding method"),
         shape = guide_legend("Model", override.aes = list(shape=c(1,2), color="black"))) +
  labs(x = NULL, y = "Classification success", title = "xgboost and rf results, midwest survey data, 10 cv folds") +
  theme_bw() + scale_x_discrete(limits = c("CA0", "CAP", "PCO"))
