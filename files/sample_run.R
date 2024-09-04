### Source attribution of the sample data ###

## Load libraries
library(tidyverse)
library(ranger)

## Load functions
source("methods/hamming.R")
source("methods/ca_unbiased.R")             # CA method with new levels scored as zero
source("methods/pco.R")                     # PCO method
source("methods/cap.R")                     # CAP method

## Load sample data 
# data is genetic data with alleles (levels) from 10 genes (variables) across 3 sources (classes), n=10 for each source
load("data/CapItOff.RData") # sample data
load("data/SeqDat.RData") # sequence information for sample data

# training data is the source data
Dat.train <- CapItOff |> filter(class != "Human") |> droplevels() |> mutate(across(everything(), factor)) 

# test data is the human data to be attributed to a source
Dat.test <- CapItOff |> filter(class == "Human") |> droplevels() |> mutate(across(everything(), factor)) 

## Prepare data
# calculate hamming distances between levels of predictor variable
genes <- names(SeqDat |> select(starts_with("Var")))
list_of_distance_matrices <- map(genes, ~ dfun(gene=., dat=CapItOff, seqdat=SeqDat))
names(list_of_distance_matrices) <- genes

# set seed
set.seed(123)

# prepare data using method of choice:
#' reminder of terminology: k = number ca axes (default is num.classes-1), 
#'                          m = number pco axes (default is num.classes-1), 
#'                          mp = propG (default is 100%), 
#'                          c = number cap axes (default is num.classes-1)
#'                          ck = number ca axes within cap (default is num.classes-1)
#'                          cm = number pco axes within cap (default is num.classes-1)
#'                          cmp = propG within cap (default is 100%)

# 1. ca unbiased method
train <- prepare_training_ca0(Dat.train, starts_with("Var"), "class")
test <- prepare_test_ca0(Dat.test, train$extra, "id")

# 2. pco method
train <- prepare_training_pco(Dat.train, starts_with("Var"), "class", d=list_of_distance_matrices, mp=99)
test <- prepare_test_pco(Dat.test, train$extra, "id")

# 3. cap method
train <- prepare_training_cap(Dat.train, starts_with("Var"), "class", d=list_of_distance_matrices, cmp=95)
test <- prepare_test_cap(Dat.test, train$extra, "id")


## generate random forest models
rf_mod <- ranger(class ~ ., data=train$training, oob.error = TRUE, num.trees=500, respect.unordered.factors = TRUE)

## make predictions for Human data
Prediction <- predict(rf_mod, data=test, predict.all = FALSE)$predictions
table(Prediction) |> as.data.frame() # counts of predictions for each source

