## Extend cerda to CAP

# Load libraries
library(here)
library(tidyverse)
library(ranger)
library(caret)   # for createFolds
library(stringdist)   # for distance matrix

# Load in our functions
source("methods/tree_predictions.R") # to pull out individual tree predictions
source("methods/ranger_mods.R")
source("methods/libs_fns.R")

# Load data
midwest_raw <- read.csv(here("cerda", "midwest","midwest_data","midwest_survey.csv"), header=T, na.strings=c(""," ","NA"))

# Data preprocessing
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

d <- list(categorical_variable = stringdistmatrix(midwest$categorical |> unique(), midwest$categorical |> unique(), 
                                                  method = "lv", useNames = "strings"))

# Make a list for the methods
list.of.data <- list(ca=midwest,
                     binary=midwest,
                     ca0=midwest,
                     pco=midwest,
                     similarity=midwest,
                     cap=midwest)


# Split data into training and test sets using folds
set.seed(3)
flds <- createFolds(y=midwest$census_region, k=10) 
list.of.train <- map(list.of.data, function(x) {map(flds, ~slice(x, {-.}))})
list.of.test <- map(list.of.data, function(x) {map(flds, ~slice(x, .))})

# Table of counts in each target region shows that the census regions are very imbalanced.
list.of.test$ca |> bind_rows() |> pull(census_region) |> table()

# Transform the predictor variable(s) according to each method
# Note that the ca method doesn't allow for numeric columns or else we could just use the ca method on the transformed data.
list.of.methods <- list(ca="ca", binary="binary", ca0="ca0", pco="pco", similarity="similarity", cap="cap")

#' reminder of terminology: k = number ca axes (default is num.classes-1), 
#'                          m = number pco axes (default is num.classes-1), 
#'                          mp = propG (default is 100%), 
#'                          c = number cap axes (default is num.classes-1)
list.of.prepped.data <- pmap(list(list.of.train, list.of.test, list.of.methods), 
                             function(x, y, z) {
                               map2(x, y, ~prep_data(.x, .y, method=z, 
                                                     class = "census_region", id="id", var_id="categorical", d=d, mp=99))})

# Train random forest models
#First merge the prepped data with the ignored (numeric) columns, then train the ranger models
numtrees <- 500
list.of.models <- pmap(
  list(list.of.prepped.data, list.of.train, list.of.test, list.of.methods), 
  function(prepped.folds, train.folds, test.folds, method){
    pmap(list(prepped.folds, train.folds, test.folds), 
         function(prepped.dat, train.dat, test.dat, ntrees=numtrees, class="census_region", id="id"){
           other.train <- train.dat |> select(-c(id, categorical_variable, census_region))
           other.test <- test.dat |> select(-c(categorical_variable, census_region))
           train_dat <- prepped.dat$train$training |> bind_cols(other.train) # hoping this lines up correctly
           test_dat <- prepped.dat$test |> left_join(other.test, by="id")
           classes <- train_dat |> pull(class) #|> droplevels()
           ranger_mod <- if(method=="binary") {
             ranger(dependent.variable.name = class, data = train_dat, 
                    oob.error = TRUE, num.trees=ntrees, respect.unordered.factors = TRUE)
           } else {
             ranger(classes ~ ., data=train_dat |> select(-any_of(class)), 
                    oob.error = TRUE, num.trees=ntrees, respect.unordered.factors = TRUE)
           }
           out <- list(extras=prepped.dat$train$extra, train=train_dat, test=test_dat, 
                       ranger_mod=ranger_mod, class=class, id=id)
           out
         })
  })


# Assess the models
#Predict the test data and merge with true class
list.of.forest_predictions <- map2(list.of.models, list.of.test, function(model, dat){  
  mod.list <- map(model, "ranger_mod")  
  prep.test.list <- map(model, "test")  
  class.list <- map(dat, function(x) x |> select("id", "census_region"))
  predictions.list <- pmap(list(mod.list, prep.test.list, class.list), 
                           function(mod, prep.test, class) {
                             predictions <- predict(mod, prep.test)$predictions  #note predict.all = FALSE 
                             bind_cols(class, prediction=predictions)}
  )}
)

#Pull out individual tree decisions for each observation, and identify observations with unique levels in the testing data.
doParallel::registerDoParallel()
list.of.tree_predictions <- map2(list.of.models, list.of.test, function(model, dat){
  mod.list <- map(model, "ranger_mod")
  prep.test.list <- map(model, "test")
  extras.list <- map(model, "extras")
  class.list <- map(model, "class")
  id.list <- map(model, "id")
  test.list <- dat # original, untransformed data
  predictions.list <- pmap(list(mod.list, prep.test.list, extras.list, class.list, id.list, test.list),
                           function(mod, prep.test, extras, class, id, test) {
                             uniques <- is_unique(data=test, list_of_extras=extras)
                             tree_preds <- predict_by_tree(mod=mod, new_data=prep.test, new_unique=uniques, id=id)
                             tree_preds |> left_join(test |> select(id, any_of(class)), by="id")}
  )}
)

# Forest misclassification rates
#Calculate the proportion of correct forest classifications (ie misclassification rate)
list.of.forest_mc <- map(list.of.forest_predictions, function(x) {
  # merge the folds
  preds <- bind_rows(x, .id="fold") 
  # weights for each fold
  w <- preds |> group_by(fold) |> tally() |> pull(n)
  # misclassification rates of the folds
  lm <- preds |> mutate(wrong = census_region != prediction) |>
    group_by(fold) |> summarise(mc = sum(wrong)/n()) |> 
    lm(mc ~ 1, weights = w, data=_) |> summary()
  # weighted mean and standard error of the misclassification rates
  av <- coef(lm)[,"Estimate"]
  se <- coef(lm)[,"Std. Error"]
  # confusion matrix
  conf <- preds |> group_by(fold, census_region, prediction) |> 
    summarise(n=n()) |> mutate(N=sum(n), p=n/sum(n)) |> ungroup() |>
    mutate(w = as.numeric(paste(factor(fold, labels=w))), 
           t_p = paste(census_region, prediction, sep="_")) |> 
    split(~t_p)
  # weighted mean and standard error of the confusion matrices
  conf.lm <- map(conf, function(x) {lm(p~1, weights=w, data=x) |> summary()})
  conf.av <- map(conf.lm, function(x) coef(x)[,"Estimate"])
  conf.se <- map(conf.lm, function(x) coef(x)[,"Std. Error"])
  out <- list(av=av, se=se, conf=conf, conf.av=conf.av, conf.se=conf.se)
  out})

# Tree misclassification rates
#Calculate misclassification rates of the individual trees according to the use of absent levels 
list.of.tree_mc <- map(list.of.tree_predictions, function(x){
  # merge the folds
  preds <- bind_rows(x, .id="fold") 
  # misclassification rates of the folds
  preds |> mutate(uses_unique = if_else(uses_unique == 0, "No", "Yes")) |>
    group_by(fold, uses_unique, census_region, prediction) |>
    summarise(n=n()) |>
    group_by(fold, uses_unique, census_region) |>
    mutate(N = sum(n)) |>
    mutate(prop = n/N)
})

# merge results
results <- map2(list.of.tree_predictions, list.of.forest_predictions, function(tree, forest){
  #tree: misclassification rate
  tree.mc <- tree |> bind_rows() |> mutate(uses_unique = if_else(uses_unique == 0, "No", "Yes")) |>
    group_by(uses_unique, census_region, prediction) |>
    summarise(n=n()) |>
    group_by(uses_unique, census_region) |>
    mutate(N = sum(n)) |>
    mutate(prop = n/N)
  tree.yes <- tree.mc |> ungroup() |> filter(uses_unique=="Yes") |> arrange(desc(prop)) |> 
    select(census_region, prediction, prop) |> 
    mutate(census_region=make.names(census_region), prediction=make.names(prediction))
  tree.no <- tree.mc |> ungroup() |> filter(uses_unique=="No") |> arrange(desc(prop)) |> 
    select(census_region, prediction, prop) |> 
    mutate(census_region=make.names(census_region), prediction=make.names(prediction))
  tree.mc <- right_join(tree.yes, tree.no, by=c("census_region", "prediction"), suffix=c(".yes",".no"))
  
  #forest: misclassification rate
  forest.mc <- forest |> bind_rows() |> group_by(census_region, prediction) |> 
    summarise(n=n()) |> mutate(N = sum(n)) |> mutate(prop = n/N) |> select(-n,-N) |> 
    mutate(census_region = make.names(census_region), prediction = make.names(prediction))
  
  #merge
  full_join(tree.mc, forest.mc, by=c("census_region", "prediction")) |> 
    replace_na(list(prop=0, prop.yes=0, prop.no=0)) |> 
    rename(tree.absent=prop.yes,tree.no.absent=prop.no,forest=prop)
})


#save(results, file="cerda/midwest/midwest_results/Cerda_cap_results_500trees_8axes.R")
load("cerda/midwest/midwest_results/Cerda_cap_results_500trees_8axes.R")  #results

# Plot results
#forest results
p1 <- results |> bind_rows(.id="method") |> 
  pivot_longer(cols=c(tree.absent, tree.no.absent, forest), names_to = "result", values_to = "p") |> 
  filter(result == "forest") |> 
  filter(prediction == census_region) |> 
  mutate(method = factor(method, levels=c("ca", "binary", "ca0", "pco", "similarity", "cap"))) |> 
  mutate(across(where(is.character), factor)) |> 
  mutate(pch = ifelse(census_region==prediction, 21,1)) |> 
  ggplot(aes(x=prediction, y=p)) + 
  geom_point(aes(colour=census_region, shape=I(pch), bg=census_region), size=3) +
  theme_bw(base_size = 11) + 
  scale_y_continuous("Proportion of predictions\n",expand=c(0.03,0.03),limits=c(0,1)) +
  #theme(legend.position = "bottom") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(strip.text.y = element_text(size = 8)) +
  facet_wrap(~method, nrow=1) +
  labs(x="\nPredicted census region", colour="True census region", bg="True census region")

p1

#tree results
p2 <- results |> bind_rows(.id="method") |> 
  pivot_longer(cols=c(tree.absent, tree.no.absent, forest), names_to = "result", values_to = "p") |> 
  filter(result != "forest") |> 
  mutate(method = factor(method, levels=c("ca", "binary", "ca0", "pco", "similarity", "cap"))) |> 
  mutate(result = factor(result, levels=c("tree.absent", "tree.no.absent"))) |> 
  mutate(across(where(is.character), factor)) |> 
  mutate(pch = ifelse(census_region==prediction, 24,1)) |> 
  ggplot(aes(x=prediction, y=p, group=result)) + 
  geom_point(aes(colour=result, shape=I(pch), bg=result), size=3) +
  theme_bw(base_size = 11) + 
  scale_y_continuous("Proportion of predictions\n",expand=c(0.03,0.03),limits=c(0,1)) +
  theme(legend.position = "bottom") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.text.y = element_text(size = 8)) +
  facet_grid(census_region~method) +
  scale_colour_manual(values=c("#9d1001","#019d10"), labels=c("Yes", "No")) + 
  scale_fill_manual(values=c("#9d1001","#019d10"), labels=c("Yes", "No")) + 
  labs(x="\nPredicted census region", colour="Absent levels used in prediction", bg="Absent levels used in prediction")

p2

png("plot1.png", width=960, height = 1024)
p1
dev.off()

png("plot2.png", width=960, height = 1024)
p2
dev.off()

