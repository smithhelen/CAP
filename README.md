CAP
================

-   <a href="#overview" id="toc-overview">Overview</a>
-   <a href="#contents-of-this-repository"
    id="toc-contents-of-this-repository">Contents of this repository</a>
-   <a href="#data-for-exploration" id="toc-data-for-exploration">Data for
    exploration</a>

### Overview

CAP combines the PCO and CA methods for random forest models

This repository is to accompany the manuscript 
**“CAP”**
by *HL Smith, PJ Biggs, NP French, ANH Smith,* and *JC Marshall* (2024).
**CAP** details new methods for encoding categorical variables for random forest predictive models which are unbiased in the presence of absent levels.

This repository contains all the code for the methods described in the paper as well as the simulation studies for reproducibility.

### Contents of this repository

The key contents are organised as follows:

-   CAP
    -   README.Rmd
    -   CAP.Rproj
    -   data
        -   NotSure.Rdata
        -   SeqDat.Rdata
    -   files
        -   sample_run.R
        -   simulation_cap
    -   methods
        -   ca.R
        -   ca_unbiased.R
        -   pco.R
        -   cap.R
        -   recipe_ca.R
        -   recipe_ca0.R
        -   recipe_pco.R
        -   recipe_cap.R
        -   hamming.R
        -   helpers.R
        -   libs_fns.R
        -   tree_predictions.R

The R project is called `CAP.Rproj`.

The directory `methods` is where all of the functions for this project
are stored.

-   `ca.R` contains the code for the original method where observations
    with absent levels are always directed to the right branch at a
    tree split
-   `ca_unbiased.R` contains the code for the adapted ca method where
    observations with absent levels are directed according to the *a
    priori* hypothesis of equal class distribution
-   `pco.R` contains the code for our new method where observations with
    absent levels are directed through the random forest according to
    their *similarity* to other levels
-   `cap.R` contains the code for our new method where observations with
    absent levels are directed through the random forest according to
    their *similarity* to other levels under the constraint that the
    direction of similarity must be in the direction of class difference
-   `hamming.R` for generating a matrix of Hamming distances between
    levels of a factor
-   `helpers.R` contains code for functions used internally in the other
    methods
-   `libs_fns.R` for easy loading of required libraries and methods
-   `tree_predictions.R` contains the code for pulling out individual
    tree predictions from a random forest
-   `recipe_ca.R` contains the ca method in tidymodels form
-   `recipe_ca0.R` contains the ca unbiased method in tidymodels form
-   `recipe_pco.R` contains the pco method in tidymodels form        
-   `recipe_cap.R` contains the cap method in tidymodels form        

The directory `files` contains the file `sample_run` which includes an
example of running the methods and generating source attribution
estimates, including individual tree predictions.

The directory `data` contains the data used in the example file
`sample_run`.

### Data for exploration

The data to accompany the example file is a small dataset for the
purpose of demonstrating the methods in the paper **CAP**.

There are two files. The file `SeqDat` contains the aligned nucleotide
sequencing data for 10 genes from each of 40 *Camplyobacter* isolates
collected from four sources (Sheep, Cattle, Chicken, and Human). The
file `BlahBlah` contains the corresponding allelic information
for the aligned sequences.
