CAP
================

-   <a href="#overview" id="toc-overview">Overview</a>
-   <a href="#contents-of-this-repository"
    id="toc-contents-of-this-repository">Contents of this repository</a>
-   <a href="#data-for-exploration" id="toc-data-for-exploration">Data for
    exploration</a>

### Overview


This repository is to accompany the manuscript 
**“To CAP it Off - Encoding Categorical Variables using Canonical Analysis of Principal Coordinates”**
by *HL Smith, PJ Biggs, NP French, ANH Smith,* and *JC Marshall* (2024).

**To CAP it Off** details new methods for encoding categorical variables for random forest predictive models which are unbiased in the presence of absent levels.
The CAP-encoding method combines the PCO-encoding and the CA-unbiased-encoding methods for encoding categorical predictor variables, including absent levels.

This repository contains all the code for the methods and simulation studies described in the paper, as well as example data to apply the methods to the source attribution of *Campylobacter* species.

### Contents of this repository

The key contents are organised as follows:

-   CAP
    -   README.Rmd
    -   CAP.Rproj
    -   data
        -   CapItOff.Rdata
        -   SeqDat.Rdata
    -   files
        -   sample_run.R
        -   to_cap_it_off_code.R
    -   methods
        -   ca_unbiased.R
        -   pco.R
        -   cap.R
        -   hamming.R
        -   helpers.R

The R project is called `CAP.Rproj`.

The directory `methods` is where all of the functions for this project
are stored.

-   `ca_unbiased.R` contains the code for the adapted ca method where
    observations with absent levels are directed according to the *a
    priori* hypothesis of equal class distribution
-   `pco.R` contains the code for our new method where observations with
    absent levels are directed through the random forest according to
    their *similarity* to other levels
-   `cap.R` contains the code for our new method where observations with
    absent levels are directed through the random forest according to
    their *similarity* to other levels under the constraint that the
    direction of similarity must be in the direction of greatest class
    difference
-   `hamming.R` for generating a matrix of Hamming distances between
    levels of a factor
-   `helpers.R` contains code for functions used internally in the other
    methods
-   `recipe_ca0.R` contains the ca unbiased method in tidymodels form
-   `recipe_pco.R` contains the pco method in tidymodels form        
-   `recipe_cap.R` contains the cap method in tidymodels form        

The directory `files` contains two files 
-   `sample_run.R` which includes an example of running the methods and generating source attribution estimates
-   `to_cap_it_off_code.R` which includes the code for the simulation studies run in the manuscript, including data generation

The directory `data` contains the two data files used in `sample_run`.

### Data for exploration

The data to accompany the example file is a small dataset for the
purpose of demonstrating the methods in the paper 
**“To CAP it Off - Encoding Categorical Variables using Canonical Analysis of Principal Coordinates”**.

There are two files:
1. `SeqDat.R` contains the aligned nucleotide
sequencing data for 10 genes from each of 40 *Camplyobacter* isolates
collected from four sources (Sheep, Cattle, Chicken, and Human)
2. `CapItOff.R` contains the corresponding allelic information
for the aligned sequences.
