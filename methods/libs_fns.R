## Load libraries and functions all in one - CAP methods ##

# Load libraries
library(tidyverse)
library(ranger)

# Load in our functions
source("methods/ca.R")                # CA method
source("methods/binary.R")            # Binary method
source("methods/ca_unbiased.R")       # CA method with new levels scored as zero
source("methods/pco.R")               # PCO method
source("methods/cap.R")               # CAP method
source("methods/similarity.R")        # similarity method
