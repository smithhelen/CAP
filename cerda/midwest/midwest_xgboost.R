## Extend cerda to xgboost

# Load libraries
library(xgboost)
library(here)
library(tidyverse)
library(ranger)
library(caret)   # for createFolds
library(stringdist)   # for distance matrix
library(tidymodels)
library(finetune)
library(vip)

# Load in our functions
source("methods/tree_predictions.R") # to pull out individual tree predictions
source("methods/ranger_mods.R")
source("methods/libs_fns.R")
source("methods/ngram.R")

# Load data
midwest_raw <- read.csv(here("cerda", "midwest","midwest_data","midwest_survey.csv"), header=T, na.strings=c(""," ","NA"))

# Data preprocessing#######################################################################################################################################################
# -- remove rows with missing values for the target variable or in any explanatory variable other than the selected categorical variable (reduces to 2421 observations);  
# -- for the selected categorical variable replace missing entries by the string ‘nan’;  
# -- transform all entries for the categorical variable to lower case;  
# -- convert all variables to ordinal, binary, or nominal factors.
# Note that Cerda 2018 standardized every column of the feature matrix to a unit variance

midwest_tidy <- midwest_raw |> 
  # remove id (because it contains rows that will be removed) and respondent
  select(-id, -respondent) |> 
  # fill in binary zeros
  mutate(across(IL:WY, \(.) {replace(.,is.na(.),0)})) |> 
  # remove rows with missing values for the target variable census_region
  filter(!is.na(census_region)) |> 
  # remove na from categorical_variable to avoid na in distance matrix
  # need to change to string (rather than "") so can use tidymodels with probability forests
  mutate(across(categorical_variable, \(.) {replace(.,is.na(.),"nan")})) |> 
  # remove rows with missing values for explanatory variables other than the selected categorical variable
  drop_na(where(is.character)) |> 
  # transform all entries for the categorical variable to lower case
  mutate(categorical_variable = tolower(categorical_variable)) |> 
  # re-add id
  rownames_to_column("id") |> 
  # convert ordinal variables to ordered factors
  mutate(education = factor(education, levels=c("Some college", "Less than high school degree", "High school degree", "Associate or bachelor degree", "Graduate degree"), ordered=TRUE, exclude = NULL), 
         age = factor(age, levels=c("18-29", "30-44", "45-60", "> 60"), ordered=TRUE, exclude = NULL), 
         household_income = factor(household_income, levels=c("$0 - $24,999", "$25,000 - $49,999", "$50,000 - $99,999", "$100,000 - $149,999", "$150,000+"), ordered=TRUE, exclude = NULL),
         personal_id = factor(personal_id, levels=c("Not at all","Not much","Some","A lot"), ordered=TRUE, exclude = NULL)) |> 
  # convert character variables to factors
  mutate(across(!where(is.numeric), factor)) |> 
  # convert all predictor variables other than the categorical_variable to numeric
  mutate(across(where(is.factor) & -c(id, census_region, categorical_variable), as.numeric))

#' Because of all the spaces and dots in the categorical variable levels, 
#' `make.names()` will reduce the number of unique levels, so we need to make them unique to start with.

unique_levels <- data.frame(categorical_variable = levels(midwest_tidy$categorical_variable)) |> 
  mutate(unique_categorical_variable = categorical_variable |> make.names(unique = TRUE))
midwest <- midwest_tidy |> left_join(unique_levels, by="categorical_variable") |> mutate(categorical_variable = as.factor(unique_categorical_variable), .keep="unused") 
str(midwest)

#The categorical variables are either binary or ordinal apart from 'categorical_variable'.
#The different methods will use different encodings of the categorical_variable column.  

# Calculate distance matrix
# Levenshtein distance "lv" ("hamming" only works when strings are the same length)
# It needs to be in a list for the pco and similarity functions though.

# 3 gram needs to have minimum of 3 characters, can't add a space as it gets removed with make.names() so add a "_"
midwest.3g <- midwest |> mutate(categorical_variable = case_when(categorical_variable == "nw" ~ "nw_", 
                                        categorical_variable == "ca" ~ "ca_", 
                                        categorical_variable == "se" ~ "se_", 
                                        categorical_variable == "ne" ~ "ne_", 
                                        categorical_variable == "sf" ~ "sf_", 
                                        categorical_variable == "wa" ~ "wa_", 
                                        categorical_variable == "ca" ~ "ca_", 
                                        TRUE ~ categorical_variable) |> as.factor())

d.lv <- list(categorical_variable = stringdistmatrix(midwest$categorical_variable |> unique(), midwest$categorical_variable |> unique(), method = "lv", useNames = "strings"))
d.3gram <- list(categorical_variable = ngram_fun("categorical_variable", midwest.3g, midwest.3g, ngram=3))
sim.3gram.cerda <- list(categorical_variable = ngram_fun.cerda("categorical_variable", midwest.3g, midwest.3g, ngram=3))


# score each method#######################################################################################################################################################
score_fn <- function(x, method, d, m=NULL, mp=99){
  #' id the training vs test set
  df = x$data
  train.df = analysis(x)
  test.df = assessment(x)
  
  # pull out categorical column and encode
  prepped_var <- prep_data(train.df, test.df, method=method, d=d, id="id", class="census_region", var_id="categorical", m=m, mp=mp)
  
  # merge with numeric columns
  train_dat <- prepped_var$train$training |> bind_cols(train.df |> select(id)) |> left_join(train.df) |> select(-categorical_variable)
  test_dat <- prepped_var$test |> left_join(test.df) |> select(-categorical_variable)
  prepped.dat <- left_join(df |> select(id), rbind(train_dat, test_dat)) # to get in original order (id is a factor)

  # return dataframe for fold
  prepped.dat
}

# need to score Cerda at start so column names are the same to use in the recipe step
cerda_data <- (prepare_training_similarity(midwest, starts_with("categorical"), class="census_region", d=d.lv))$training |> bind_cols(midwest |> select(id)) |> left_join(midwest) |> select(-categorical_variable)
cerda_3g_data <- (prepare_training_similarity(midwest.3g, starts_with("categorical"), class="census_region", d=sim.3gram.cerda))$training |> bind_cols(midwest.3g |> select(id)) |> left_join(midwest.3g) |> select(-categorical_variable)

#separate into folds########################################################################################################################

# separate into folds
set.seed(246)
folds <- vfold_cv(midwest, strata = census_region)
folds.3g <- vfold_cv(midwest.3g, strata = census_region)
folds.cerda <- vfold_cv(cerda_data, strata = census_region)
folds.cerda.3g <- vfold_cv(cerda_3g_data, strata = census_region)



## CHECK
#x <- folds$splits[[1]]
#method = "ca0"
#prepped.dat <- score_fn(x, method, d=NULL)
#folds$splits[[1]]$data <- prepped.dat
#analysis(folds$splits[[1]])
#assessment(folds$splits[[1]])

# replace folds
ca0_folds <- folds
for(i in 1:10){
  x <- ca0_folds$splits[[i]]
  prepped.dat <- score_fn(x, method="ca0", d=NULL)
  ca0_folds$splits[[i]]$data <- prepped.dat
}

cerda_folds <- folds.cerda  # have to prescore everything

pco_folds <- folds
for(i in 1:10){
  x <- pco_folds$splits[[i]]
  prepped.dat <- score_fn(x, method="pco", d=d.lv, mp=95)
  pco_folds$splits[[i]]$data <- prepped.dat
}

# pull out number of dimensions that was the largest
m=pco_folds$splits |> map(function(x) x$data |> ncol()) |> unlist() |> max() - (midwest |> ncol() -1)
pco_folds <- folds
for(i in 1:10){
  x <- pco_folds$splits[[i]]
  prepped.dat <- score_fn(x, method="pco", d=d.lv, m=m, mp=NULL)
  pco_folds$splits[[i]]$data <- prepped.dat
}

cap_folds <- folds
for(i in 1:10){
  x <- cap_folds$splits[[i]]
  prepped.dat <- score_fn(x, method="cap", m=m, mp=NULL, d=d.lv)
  cap_folds$splits[[i]]$data <- prepped.dat
}

cerda_3gram_folds <- folds.cerda.3g  # have to prescore everything

pco_3gram_folds <- folds.3g
for(i in 1:10){
  x <- pco_3gram_folds$splits[[i]]
  prepped.dat <- score_fn(x, method="pco", d=d.3gram, mp=95)
  pco_3gram_folds$splits[[i]]$data <- prepped.dat
}

# pull out number of dimensions that was the largest
m3=pco_3gram_folds$splits |> map(function(x) x$data |> ncol()) |> unlist() |> max() - (midwest |> ncol() -1)
pco_3gram_folds <- folds.3g
for(i in 1:10){
  x <- pco_3gram_folds$splits[[i]]
  prepped.dat <- score_fn(x, method="pco", d=d.3gram, m=m3, mp=NULL)
  pco_3gram_folds$splits[[i]]$data <- prepped.dat
}

cap_3gram_folds <- folds.3g
for(i in 1:10){
  x <- cap_3gram_folds$splits[[i]]
  prepped.dat <- score_fn(x, method="cap", d=d.3gram, m=m3, mp=NULL)
  cap_3gram_folds$splits[[i]]$data <- prepped.dat
}

#specify models#######################################################################################################################################################

# specify models
xg_spec <- 
  boost_tree(trees = tune(), #default 15
             min_n = tune(), #default 1
             mtry = tune(), #default 
             #tree_depth = tune(), #default 6           
             learn_rate = 0.01) |>
  set_engine("xgboost") |> 
  set_mode("classification")

xg_recipe.ca0 <- recipe(census_region ~ ., data = ca0_folds$splits[[1]]$data) |> 
  update_role(id, new_role = "id")
xg_workflow.ca0 <- workflow(xg_recipe.ca0, xg_spec)

xg_recipe.cerda <- recipe(census_region ~ ., data = cerda_data) |> 
  update_role(id, new_role = "id")
xg_workflow.cerda <- workflow(xg_recipe.cerda, xg_spec)

xg_recipe.pco <- recipe(census_region ~ ., data = pco_folds$splits[[1]]$data) |> 
  update_role(id, new_role = "id")
xg_workflow.pco <- workflow(xg_recipe.pco, xg_spec)

xg_recipe.cap <- recipe(census_region ~ ., data = cap_folds$splits[[1]]$data) |> 
  update_role(id, new_role = "id")
xg_workflow.cap <- workflow(xg_recipe.cap, xg_spec)

xg_recipe.cerda.3g <- recipe(census_region ~ ., data = cerda_3g_data) |> 
  update_role(id, new_role = "id")
xg_workflow.cerda.3g <- workflow(xg_recipe.cerda.3g, xg_spec)

xg_recipe.pco.3g <- recipe(census_region ~ ., data = pco_3gram_folds$splits[[1]]$data) |> 
  update_role(id, new_role = "id")
xg_workflow.pco.3g <- workflow(xg_recipe.pco.3g, xg_spec)

xg_recipe.cap.3g <- recipe(census_region ~ ., data = cap_3gram_folds$splits[[1]]$data) |> 
  update_role(id, new_role = "id")
xg_workflow.cap.3g <- workflow(xg_recipe.cap.3g, xg_spec)

#Tune models########################################################################################################################

## Tune models
doParallel::registerDoParallel()
set.seed(345)
xg_resamples.ca0 <- tune_race_anova(
  xg_workflow.ca0,
  resamples = ca0_folds,
  grid = 15,
  control = control_race(verbose_elim = TRUE)
)

xg_resamples.cerda <- tune_race_anova(
  xg_workflow.cerda,
  resamples = cerda_folds,
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

xg_resamples.cerda.3g <- tune_race_anova(
  xg_workflow.cerda.3g,
  resamples = cerda_3gram_folds,
  grid = 15,
  control = control_race(verbose_elim = TRUE)
)

xg_resamples.pco.3g <- tune_race_anova(
  xg_workflow.pco.3g,
  resamples = pco_3gram_folds,
  grid = 15,
  control = control_race(verbose_elim = TRUE)
)

xg_resamples.cap.3g <- tune_race_anova(
  xg_workflow.cap.3g,
  resamples = cap_3gram_folds,
  grid = 15,
  control = control_race(verbose_elim = TRUE)
)


#find best parameters ########################################################################################################################
## find best parameters
collect_metrics(xg_resamples.ca0)
plot_race(xg_resamples.ca0)
show_best(xg_resamples.ca0, "accuracy")
best_ca0 <- select_best(xg_resamples.ca0, "accuracy")
xg_final_workflow.ca0 <- finalize_workflow(xg_workflow.ca0, best_ca0)
  
collect_metrics(xg_resamples.cerda)
plot_race(xg_resamples.cerda)
show_best(xg_resamples.cerda, "accuracy")
best_cerda <- select_best(xg_resamples.cerda, "accuracy")
xg_final_workflow.cerda <- finalize_workflow(xg_workflow.cerda, best_cerda)

collect_metrics(xg_resamples.pco)
plot_race(xg_resamples.pco)
show_best(xg_resamples.pco, "accuracy")
best_pco <- select_best(xg_resamples.pco, "accuracy")
xg_final_workflow.pco <- finalize_workflow(xg_workflow.pco, best_pco)

collect_metrics(xg_resamples.cap)
plot_race(xg_resamples.cap)
show_best(xg_resamples.cap, "accuracy")
best_cap <- select_best(xg_resamples.cap, "accuracy")
xg_final_workflow.cap <- finalize_workflow(xg_workflow.cap, best_cap)

collect_metrics(xg_resamples.cerda.3g)
plot_race(xg_resamples.cerda.3g)
show_best(xg_resamples.cerda.3g, "accuracy")
best_cerda.3g <- select_best(xg_resamples.cerda.3g, "accuracy")
xg_final_workflow.cerda.3g <- finalize_workflow(xg_workflow.cerda.3g, best_cerda.3g)

collect_metrics(xg_resamples.pco.3g)
plot_race(xg_resamples.pco.3g)
show_best(xg_resamples.pco.3g, "accuracy")
best_pco.3g <- select_best(xg_resamples.pco.3g, "accuracy")
xg_final_workflow.pco.3g <- finalize_workflow(xg_workflow.pco.3g, best_pco.3g)

collect_metrics(xg_resamples.cap.3g)
plot_race(xg_resamples.cap.3g)
show_best(xg_resamples.cap.3g, "accuracy")
best_cap.3g <- select_best(xg_resamples.cap.3g, "accuracy")
xg_final_workflow.cap.3g <- finalize_workflow(xg_workflow.cap.3g, best_cap.3g)

all_workflows <- list(ca0=xg_final_workflow.ca0, 
                      cerda=xg_final_workflow.cerda, 
                      pco=xg_final_workflow.pco, 
                      cap=xg_final_workflow.cap, 
                      cerda_3g=xg_final_workflow.cerda.3g, 
                      pco_3g=xg_final_workflow.pco.3g, 
                      cap_3g=xg_final_workflow.cap.3g)
save(all_workflows, file = "final_workflows_midwest_mp95.Rdata")

#separate into folds for a second time########################################################################################################################

# separate into folds for a second time
set.seed(456)
folds <- vfold_cv(midwest, strata = census_region)
folds.3g <- vfold_cv(midwest.3g, strata = census_region)
folds.cerda <- vfold_cv(cerda_data, strata = census_region)
folds.cerda.3g <- vfold_cv(cerda_3g_data, strata = census_region)

# score each method
ca0_folds <- folds
for(i in 1:10){
  x <- ca0_folds$splits[[i]]
  prepped.dat <- score_fn(x, method="ca0", d=NULL)
  ca0_folds$splits[[i]]$data <- prepped.dat
}

cerda_folds <- folds.cerda  # have to prescore everything

pco_folds <- folds
for(i in 1:10){
  x <- pco_folds$splits[[i]]
  prepped.dat <- score_fn(x, method="pco", d=d.lv, mp=95)
  pco_folds$splits[[i]]$data <- prepped.dat
}

# pull out number of dimensions that was the largest
m=pco_folds$splits |> map(function(x) x$data |> ncol()) |> unlist() |> max() - (midwest |> ncol() -1)
pco_folds <- folds
for(i in 1:10){
  x <- pco_folds$splits[[i]]
  prepped.dat <- score_fn(x, method="pco", d=d.lv, m=m, mp=NULL)
  pco_folds$splits[[i]]$data <- prepped.dat
}

cap_folds <- folds
for(i in 1:10){
  x <- cap_folds$splits[[i]]
  prepped.dat <- score_fn(x, method="cap", m=m, mp=NULL, d=d.lv)
  cap_folds$splits[[i]]$data <- prepped.dat
}

cerda_3gram_folds <- folds.cerda.3g  # have to prescore everything

pco_3gram_folds <- folds.3g
for(i in 1:10){
  x <- pco_3gram_folds$splits[[i]]
  prepped.dat <- score_fn(x, method="pco", d=d.3gram, mp=95)
  pco_3gram_folds$splits[[i]]$data <- prepped.dat
}

# pull out number of dimensions that was the largest
m3=pco_3gram_folds$splits |> map(function(x) x$data |> ncol()) |> unlist() |> max() - (midwest |> ncol() -1)

pco_3gram_folds <- folds.3g
for(i in 1:10){
  x <- pco_3gram_folds$splits[[i]]
  prepped.dat <- score_fn(x, method="pco", d=d.3gram, mp=NULL, m=m3)
  pco_3gram_folds$splits[[i]]$data <- prepped.dat
}

cap_3gram_folds <- folds.3g
for(i in 1:10){
  x <- cap_3gram_folds$splits[[i]]
  prepped.dat <- score_fn(x, method="cap", d=d.3gram, m=m3, mp=NULL)
  cap_3gram_folds$splits[[i]]$data <- prepped.dat
}



#fit resamples#######################################################################################################################################################

# fit resamples
set.seed(567)
xg_rs.ca0 <- fit_resamples(
  xg_final_workflow.ca0,
  resamples = ca0_folds,
  control = control_resamples(save_pred = TRUE)
)
collect_metrics(xg_rs.ca0)

xg_rs.cerda <- fit_resamples(
  xg_final_workflow.cerda,
  resamples = cerda_folds,
  control = control_resamples(save_pred = TRUE)
)
collect_metrics(xg_rs.cerda)

xg_rs.pco <- fit_resamples(
  xg_final_workflow.pco,
  resamples = pco_folds,
  control = control_resamples(save_pred = TRUE)
)
collect_metrics(xg_rs.pco)

xg_rs.cap <- fit_resamples(
  xg_final_workflow.cap,
  resamples = cap_folds,
  control = control_resamples(save_pred = TRUE)
)
collect_metrics(xg_rs.cap)

xg_rs.cerda.3g <- fit_resamples(
  xg_final_workflow.cerda.3g,
  resamples = cerda_3gram_folds,
  control = control_resamples(save_pred = TRUE)
)
collect_metrics(xg_rs.cerda.3g)

xg_rs.pco.3g <- fit_resamples(
  xg_final_workflow.pco.3g,
  resamples = pco_3gram_folds,
  control = control_resamples(save_pred = TRUE)
)
collect_metrics(xg_rs.pco.3g)

xg_rs.cap.3g <- fit_resamples(
  xg_final_workflow.cap.3g,
  resamples = cap_3gram_folds,
  control = control_resamples(save_pred = TRUE)
)
collect_metrics(xg_rs.cap.3g)


all_metrics <- list(ca0=xg_rs.ca0, 
                      cerda=xg_rs.cerda, 
                      pco=xg_rs.pco, 
                      cap=xg_rs.cap, 
                      cerda_3g=xg_rs.cerda.3g, 
                      pco_3g=xg_rs.pco.3g, 
                      cap_3g=xg_rs.cap.3g)
save(all_metrics, file = "xg_metrics_midwest_95.Rdata")

all_metrics |> map_df(collect_metrics, .id = "method") |> 
  filter(.metric == "accuracy") |> 
  mutate(pch = case_when(method |> str_detect("3g") ~ 1,
                         method |> str_detect("ca0") ~ 15,
                         TRUE ~ 19), method = factor(gsub("_3g","",method))) |> 
  ggplot() +
  geom_point(aes(x=as.numeric(method), y=mean, col=method, shape=pch), size=4) +
  scale_shape_identity(guide = "legend", breaks=c(1,15,19), labels=c("3grams", "none", "Levenshtein")) +
  guides(shape = guide_legend("Distance measure"), col = guide_legend("Method")) +
  labs(x = NULL, y = "Classification success", title = "xgboost results, midwest survey data, 10 cv folds") +
  theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  


#Compare with random forest#######################################################################################################################################################
## Compare with random forest
rf_spec <- 
  rand_forest(trees = 500)  |>  
  set_engine("ranger", respect.unordered.factors = TRUE) |> 
  set_mode("classification")

rf_recipe.ca0 <- recipe(census_region ~ ., data = ca0_folds$splits[[1]]$data) |> 
  update_role(id, new_role = "id")
rf_workflow.ca0 <- workflow(rf_recipe.ca0, rf_spec)

rf_recipe.cerda <- recipe(census_region ~ ., data = cerda_data) |> 
  update_role(id, new_role = "id")
rf_workflow.cerda <- workflow(rf_recipe.cerda, rf_spec)

rf_recipe.pco <- recipe(census_region ~ ., data = pco_folds$splits[[1]]$data) |> 
  update_role(id, new_role = "id")
rf_workflow.pco <- workflow(rf_recipe.pco, rf_spec)

rf_recipe.cap <- recipe(census_region ~ ., data = cap_folds$splits[[1]]$data) |> 
  update_role(id, new_role = "id")
rf_workflow.cap <- workflow(rf_recipe.cap, rf_spec)

rf_recipe.cerda.3g <- recipe(census_region ~ ., data = cerda_3g_data) |> 
  update_role(id, new_role = "id")
rf_workflow.cerda.3g <- workflow(rf_recipe.cerda.3g, rf_spec)

rf_recipe.pco.3g <- recipe(census_region ~ ., data = pco_3gram_folds$splits[[1]]$data) |> 
  update_role(id, new_role = "id")
rf_workflow.pco.3g <- workflow(rf_recipe.pco.3g, rf_spec)

rf_recipe.cap.3g <- recipe(census_region ~ ., data = cap_3gram_folds$splits[[1]]$data) |> 
  update_role(id, new_role = "id")
rf_workflow.cap.3g <- workflow(rf_recipe.cap.3g, rf_spec)

set.seed(687)
rf_rs.ca0 <- fit_resamples(
  rf_workflow.ca0,
  resamples = ca0_folds,
  control = control_resamples(save_pred = TRUE)
)
collect_metrics(rf_rs.ca0)

rf_rs.cerda <- fit_resamples(
  rf_workflow.cerda,
  resamples = cerda_folds,
  control = control_resamples(save_pred = TRUE)
)
collect_metrics(rf_rs.cerda)

rf_rs.pco <- fit_resamples(
  rf_workflow.pco,
  resamples = pco_folds,
  control = control_resamples(save_pred = TRUE)
)
collect_metrics(rf_rs.pco)

rf_rs.cap <- fit_resamples(
  rf_workflow.cap,
  resamples = cap_folds,
  control = control_resamples(save_pred = TRUE)
)
collect_metrics(rf_rs.cap)

rf_rs.cerda.3g <- fit_resamples(
  rf_workflow.cerda.3g,
  resamples = cerda_3gram_folds,
  control = control_resamples(save_pred = TRUE)
)
collect_metrics(rf_rs.cerda.3g)

rf_rs.pco.3g <- fit_resamples(
  rf_workflow.pco.3g,
  resamples = pco_3gram_folds,
  control = control_resamples(save_pred = TRUE)
)
collect_metrics(rf_rs.pco.3g)

rf_rs.cap.3g <- fit_resamples(
  rf_workflow.cap.3g,
  resamples = cap_3gram_folds,
  control = control_resamples(save_pred = TRUE)
)
collect_metrics(rf_rs.cap.3g)


all_metrics <- list(ca0=list(rf=rf_rs.ca0, xg=xg_rs.ca0), 
                    cerda=list(rf=rf_rs.cerda, xg=xg_rs.cerda), 
                    pco=list(rf=rf_rs.pco, xg=xg_rs.pco), 
                    cap=list(rf=rf_rs.cap,xg=xg_rs.cap), 
                    cerda_3g=list(rf=rf_rs.cerda.3g,xg=xg_rs.cerda.3g), 
                    pco_3g=list(rf=rf_rs.pco.3g, xg=xg_rs.pco.3g), 
                    cap_3g=list(rf=rf_rs.cap.3g,xg=xg_rs.cap.3g)
                    )
save(all_metrics, file = "all_metrics_midwest_95.Rdata")

p <- all_metrics |> map_df(function(x) x |> map_df(collect_metrics, .id = "model"), .id = "method" ) |> 
  filter(.metric == "accuracy") |> 
  mutate(pch = case_when(model == "rf" ~ 21,
                        model == "xg" ~ 24), 
         dist = case_when(method |> str_detect("3g") ~ 1,
                           method |> str_detect("ca0") ~ 0.05,
                           TRUE ~ 0.5),
         method2 = factor(method),
         method = factor(gsub("_3g","",method))) 
p |> 
  ggplot() +
  geom_point(aes(x=as.numeric(method), y=mean, col=method, fill=interaction(method, dist), shape=pch), size=3, stroke=1.5) +
  scale_shape_identity(guide = "legend", breaks=c(21,24), labels=c("random forest","xgboost")) +
  scale_alpha_identity(guide = "legend", breaks=c(0.05,0.5,1), labels=c("none", "Levenshtein", "3grams")) +
  scale_colour_manual(values = c("#4C4C53", "#115896", "#7C6C65", "#BA2F00")) +
  scale_fill_manual(values = c("white",  "#7ea6c9", "#b5aeaa","#da7552", "#115896","#7C6C65", "#BA2F00")) +
  guides(col = guide_legend("Encoding method"),
         fill = "none",
         #fill = guide_legend("Distance measure"),
         shape = guide_legend("Model", override.aes = list(shape=c(1,2), color="black"))) +
  labs(x = NULL, y = "Classification success", title = "xgboost and rf results, midwest survey data, 10 cv folds", subtitle = "shading represents distance measure, none<lv<3gram") +
  theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  

