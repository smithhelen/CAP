# testing recipe step
library(tidyverse)
library(recipes)
library(rsample)
source('methods/recipe_ca.R')
source('methods/recipe_ca0.R')


load("../CAP_data/data/cgMLST_dat.RData") # SACNZ cgMLST data set (jejuni and coli)
Dat_jc <- cgMLST %>% filter(Source != "Human") %>% droplevels() %>% mutate(across(everything(), factor)) 

set.seed(3)

split <- initial_split(Dat_jc)
jc_train <- training(split)
jc_test  <- testing(split)

#### ca with scores
my_recipe <- 
  recipe(Source ~ ., data=jc_train) |>
  step_ca(starts_with("CAMP"), k=3)

prepped_recipe <- my_recipe |>
  prep()

# OK, now apply to our training data
baked_train <- prepped_recipe |>
  bake(jc_train)
baked_test <- prepped_recipe |>
  bake(jc_test)

# compare with what we had before...
source('methods/ca.R')

foo <- prepare_training_ca(jc_train, starts_with("CAMP"), "Source", k=3)
foo_train <- foo$training
switch_name <- function(nm) {
  if(nm %in% c("Source", "LabID")){return(nm)} else {
    root <- substring(nm, 1, 8)
    colnum <- substring(nm, 11)
    paste(root, "ca", colnum, sep="_")
  }
}
names(foo_train) <- map_chr(names(foo_train), switch_name)
baked_train |> anti_join(foo_train)
foo_train |> anti_join(baked_train)
foo_train |> as_tibble(); baked_train
# YAY, they're the same! :)

foo_test <- prepare_test_ca(jc_test, foo$extra, "LabID")
names(foo_test) <- map_chr(names(foo_test), switch_name)
baked_test |> anti_join(foo_test)
foo_test |> anti_join(baked_test)
foo_test |> as_tibble(); baked_test
# YAY, they're the same! :)


#### Now the unbiased version (try with different values of k as well)
my_recipe <- 
  recipe(Source ~ ., data=jc_train) |>
  step_ca_unbiased(starts_with("CAMP"), k=1)

prepped_recipe <- my_recipe |>
  prep()

# OK, now apply to our training data
baked_train <- prepped_recipe |>
  bake(jc_train)
baked_test <- prepped_recipe |>
  bake(jc_test)

# compare with what we had before...
source('methods/ca_unbiased.R')

foo <- prepare_training_ca0(jc_train, starts_with("CAMP"), "Source", k=1)
foo_train <- foo$training
switch_name <- function(nm) {
  if(nm %in% c("Source", "LabID")){return(nm)} else {
    root <- substring(nm, 1, 8)
    colnum <- substring(nm, 11)
    paste(root, "ca0", colnum, sep="_")
  }
}
names(foo_train) <- map_chr(names(foo_train), switch_name)
baked_train |> anti_join(foo_train)
foo_train |> anti_join(baked_train)
foo_train |> as_tibble(); baked_train
# YAY, they're the same! :)

foo_test <- prepare_test_ca0(jc_test, foo$extra, "LabID")
names(foo_test) <- map_chr(names(foo_test), switch_name)
baked_test |> anti_join(foo_test)
foo_test |> anti_join(baked_test)
foo_test |> as_tibble(); baked_test
# they're kinda the same


#### NOW PCO...
source("methods/recipe_pco.R")

load("../CAP_data/data/list_of_distance_matrices_all.RData")
# RIGHT, let's mess with one of them to have variance zero...
jc_train$CAMP0001 <- rep(jc_train$CAMP0001[1], nrow(jc_train))
my_recipe <- 
  recipe(Source ~ ., data=jc_train) |>
  step_pco(starts_with("CAMP"), distances = list_of_distance_matrices_all, axes=10)

prepped_recipe <- my_recipe |>
  prep()

baked_train <- prepped_recipe |>
  bake(jc_train)

# Try the old method
source("methods/pco.R")

foo <- prepare_training_pco(jc_train, starts_with("CAMP"), "Source", list_of_distance_matrices_all, m=10)$training
switch_name <- function(nm) {
  root <- substring(nm, 1, 8)
  colnum <- substring(nm, 11)
  paste(root, "pco", colnum, sep="_")
}
names(foo) <- map_chr(names(foo), switch_name)

baked_train |> anti_join(foo) # SAME! :)
