## calculate weighted Hamming distance between levels of predictor variable ####

# Load libraries
library(tidyverse)
library(seqinr) # for Hamming distances
library(stringdist) # for Hamming distances

# function to create distance (Hamming) matrices between alleles for each gene ie seq differences on alleles not individuals
# scaled by their Phandango score to allow for recombination
dfun_Phandango <- function(gene, Phandango, dat, seqdat){  
  
  ## gene stuff
  # pull out allele names for the selected gene (factor)
  alleles <- dat |> pull(gene) |> unique()
  # pull out sequence data for the alleles of the selected gene (vector of strings)
  gene_seq <- seqdat |> select(all_of(gene)) |> distinct() |> pull(gene)
  # convert strings to character vectors: (list of character vectors)
  gene_seq_char <- gene_seq |> str_split("")
  # identify matching nt positions between the different alleles of the selected gene 
  # returns a list of lists of logical vectors (TRUE/FALSEs)
  gene_allele_match <- map(gene_seq_char, 
                           function(z) {out <- map(gene_seq_char, function(y) z!=y); 
                           names(out) = alleles; out}) 
  names(gene_allele_match) <- alleles
  
  # pull out the Phandango scores for the selected gene
  score <- Phandango |> pull(score)
  
  ## calculate distance matrix
  # multiply the scores by the TRUE/FALSE vector and sum
  # gives an allele by allele dataframe of the modified Hamming distances
  d <- map_df(gene_allele_match, function(x) {map(x, function(y) {(y * score) |> sum()})}) |> 
    as.data.frame()
  rownames(d) <- colnames(d)
  return(d)
}


########################################################################################################################################################
# load mapping data and scores (mapping_score)
load(file="Phandango/Phandango_data/mapping_score.Rdata") 
# we are interested in the score column (the score.s is multiplied by the number of isolates which are affected)

# load allele and sequence data
load("../CAP_data/data/cgMLST_dat.RData") # SACNZ cgMLST data set (cgMLST)
load("../CAP_data/data/cgMLST_seq_dat.RData") # SACNZ cgMLST sequence data set (cgMLST_seq_dat)

# prepare data - pull out gene names and select animal only isolates
genes <- names(cgMLST_seq_dat %>% select(starts_with("CAMP")))
dat_alleles <- cgMLST %>% filter(Source != "Human") %>% droplevels()
dat_seq <- cgMLST_seq_dat  %>% select(-LabID) %>% filter(Source != "Human") %>% droplevels() %>% as.data.frame()

# replace na with zero for mapping_score score and separate into a list according to gene
scores_list <- mapping_score |> mutate(score = replace_na(score,0)) |> 
  group_split(gene) |> 
  set_names(mapping_score$gene |> unique() |> sort())

# calculate hamming distances between levels of predictor variable
list_of_distance_matrices_Phandango <- map2(genes, scores_list, ~ dfun_Phandango(gene=.x, Phandango=.y, dat=dat_alleles, seqdat=dat_seq))
names(list_of_distance_matrices_Phandango) <- genes
list_of_distance_matrices_Phandango <- map(list_of_distance_matrices_Phandango, as.matrix)
#save(list_of_distance_matrices_Phandango, file="Phandango/Phandango_data/list_of_distance_matrices_Phandango.RData") # run 7/4/23

# now load
load("Phandango/Phandango_data/list_of_distance_matrices_Phandango.RData")

