# testing recipe step
library(tidyverse)
library(recipes)
library(rsample)
source('methods/recipe_ca.R')
source('methods/recipe_ca_unbiased.R')


load("data/cgMLST_dat.RData") # SACNZ cgMLST data set (jejuni and coli)
Dat_jc <- cgMLST %>% filter(Source != "Human") %>% droplevels() %>% mutate(across(everything(), factor)) 

set.seed(3)

split <- initial_split(Dat_jc)
jc_train <- training(split)
jc_test  <- testing(split)

my_recipe <- 
  recipe(Source ~ ., data=jc_train) |>
  step_ca_rank(starts_with("CAMP"))

prepped_recipe <- my_recipe |>
  prep()

# OK, now apply to our training data
baked_train <- prepped_recipe |>
  bake(jc_train)

# compare with what we had before...
source('methods/ca.R')

foo <- prepare_training_ca(jc_train, starts_with("CAMP"), "Source")$training

baked_train |> anti_join(foo)
# YAY, they're the same! :)

# Now the unbiased version
my_recipe <- 
  recipe(Source ~ ., data=jc_train) |>
  step_ca_unbiased(starts_with("CAMP"))

prepped_recipe <- my_recipe |>
  prep()

# OK, now apply to our training data
baked_train <- prepped_recipe |>
  bake(jc_train)

# compare with what we had before...
source('methods/ca_unbiased.R')

foo <- prepare_training_ca0(jc_train, starts_with("CAMP"), "Source")$training

baked_train |> anti_join(foo)
# YAY, they're the same! :)

foo |> anti_join(baked_train) # Slightly different. Prob because of the way we calculate ranks?


# NOW PCO...
source("methods/recipe_pco.R")

load("data/list_of_distance_matrices_all.RData")
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

foo <- prepare_training_pco(jc_train, starts_with("CAMP"), "Source", list_of_distance_matrices_all, axes=10)$training
switch_name <- function(nm) {
  root <- substring(nm, 1, 8)
  colnum <- substring(nm, 11)
  paste(root, "pco", colnum, sep="_")
}
names(foo) <- map_chr(names(foo), switch_name)

baked_train |> anti_join(foo) # SAME! :)
