## Load libraries and functions all in one - CAP methods ##

# Load libraries
library(tidyverse)
library(ranger)
library(caret)   # for createFolds
library(varhandle) # for to.dummy()
library(recipes)
library(rsample) # for initial_split()
library(tidymodels)

# Load in our functions
source("methods/ca.R")                # CA method
source("methods/ca_unbiased.R")       # CA method with new levels scored as zero
source("methods/pco.R")               # PCO method
source("methods/cap.R")               # CAP method
source("methods/similarity.R")        # similarity method
source('methods/recipe_ca.R')         # CA recipe
source('methods/recipe_ca0.R')        # CA unbiased recipe
source("methods/recipe_pco.R")        # PCO recipe
source("methods/recipe_cap.R")        # CAP recipe
