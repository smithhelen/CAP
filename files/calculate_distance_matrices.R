# Calculate distance matrices

## Load libraries and functions ##
library(tidyverse)
source("methods/hamming.R")
source("methods/ngram.R")

## load data
load("../CAP_data/data/cgMLST_dat.RData") # SACNZ cgMLST data set (cgMLST)
load("../CAP_data/data/cgMLST_seq_dat.RData") # SACNZ cgMLST sequence data set (cgMLST_seq_dat)

# select animal isolates only
dat <- cgMLST |> filter(Source != "Human") |> droplevels()
seqdat <- cgMLST_seq_dat |> select(-LabID) |> filter(Source != "Human") |> droplevels() |> as.data.frame()
genes <- names(cgMLST_seq_dat |> select(starts_with("CAMP")))


# Hamming distance
list_of_distance_matrices <- map(genes, function(x) dfun(x, dat, seqdat))
names(list_of_distance_matrices) <- genes
#save(list_of_distance_matrices, file="../CAP_Data/data/list_of_distance_matrices.RData") # distance matrices of hamming distances among unique alleles for each gene

# ngrams
list_of_6grams <- map(genes, function(x) ngram_fun(x, dat, seqdat, ngram=6))
names(list_of_6grams) <- genes
list_of_8grams <- map(genes, function(x) ngram_fun(x, dat, seqdat, ngram=8))
names(list_of_8grams) <- genes
list_of_10grams <- map(genes, function(x) ngram_fun(x, dat, seqdat, ngram=10))
names(list_of_10grams) <- genes
list_of_12grams <- map(genes, function(x) ngram_fun(x, dat, seqdat, ngram=12))
names(list_of_12grams) <- genes
#save(list_of_6grams, file="../CAP_Data/data/list_of_6grams.RData") 
#save(list_of_8grams, file="../CAP_Data/data/list_of_8grams.RData") 
#save(list_of_10grams, file="../CAP_Data/data/list_of_10grams.RData") 
#save(list_of_12grams, file="../CAP_Data/data/list_of_12grams.RData") 



# confirm the number of ngrams to keep
list_of_3grams <- map(genes, function(x) ngram_fun(x, dat, seqdat, ngram=3))
names(list_of_3grams) <- genes
list_of_4grams <- map(genes, function(x) ngram_fun(x, dat, seqdat, ngram=4))
names(list_of_4grams) <- genes
list_of_5grams <- map(genes, function(x) ngram_fun(x, dat, seqdat, ngram=5))
names(list_of_5grams) <- genes
list_of_7grams <- map(genes, function(x) ngram_fun(x, dat, seqdat, ngram=7))
names(list_of_7grams) <- genes
list_of_9grams <- map(genes, function(x) ngram_fun(x, dat, seqdat, ngram=9))
names(list_of_9grams) <- genes
list_of_11grams <- map(genes, function(x) ngram_fun(x, dat, seqdat, ngram=11))
names(list_of_11grams) <- genes

check <- list("3grams"=list_of_3grams,
              "4grams"=list_of_4grams,
              "5grams"=list_of_5grams,
              "6grams"=list_of_6grams,
              "7grams"=list_of_7grams,
              "8grams"=list_of_8grams,
              "9grams"=list_of_9grams,
              "10grams"=list_of_10grams,
              "11grams"=list_of_11grams,
              "12grams"=list_of_12grams)
MLST <- c("CAMP0075", "CAMP0093", "CAMP0367", "CAMP0399", "CAMP0645", "CAMP1541", "CAMP1576")

mlst_genes <- map(MLST, function(x) map(check, x))
names(mlst_genes) <- MLST

mlst_dist <- map(mlst_genes, function(gene) {
  df <- map(gene, function(ngram) {
    ngram |> dist() |> as.data.frame()
    }) |> reduce(bind_cols)
  colnames(df) <- names(gene)
  df})

p <- do.call(rbind, mlst_dist) |> rownames_to_column("gene") |> mutate(gene = str_sub(gene,1,8)) |> 
  pivot_longer(cols = -gene, names_to = "ngram", values_to = "dist")

p |> mutate(ngram = factor(ngram, levels=c("3grams","4grams","5grams","6grams","7grams","8grams","9grams","10grams","11grams","12grams"))) |> ggplot(aes(y=dist, colour=ngram)) + geom_boxplot() + facet_wrap(~gene) + theme_bw()
