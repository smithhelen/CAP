### Convert NT data to aa data ###

## 1. load libraries ##
library(tidyverse)
library(seqinr)  # for s2c and getTrans
library(stringdist) # for Hamming distances


## 2. define functions ##
# function to substring the sequence into codon triplets
get_codon <- function(x) {
  substring(x,seq(1,nchar(x),3),seq(3,nchar(x)+3-1,3))
}

# function to calculate Hamming distance between aa 'alleles'
dfun <- function(gene){
  alleles <- gene |> arrange(allele) |> pull(allele) |> unique()
  x <- gene |> arrange(allele) |> distinct(aa) |> pull(aa)
  d <- stringdistmatrix(x, x, method="hamming")
  rownames(d) <- alleles
  colnames(d) <- alleles
  return(d)
}


## 3. load data ##
# sequence data
load("../CAP_data/data/cgMLST_seq_dat.RData ")


## 4. convert to amino acids and aa 'alleles' ##
# find amino acids
aa_dat <- cgMLST_seq_dat |> 
  pivot_longer(cols=starts_with("CAMP"), names_to = "gene", values_to = "sequence") |> 
#  slice(c(1:3,1344:1346)) |> # test on a few rows
  mutate(NT = map(sequence, s2c)) |> # convert sequence string into vector of characters
  mutate(codon = map(sequence, get_codon)) |> 
  mutate(aa = map(NT, getTrans)) |> 
  select(LabID, Source, gene, aa) 

# save file
#save(aa_dat, file="../CAP_data/data/aa_dat.RData")
#load("../CAP_data/data/aa_dat.RData")

# convert to alleles
allele_dat <- aa_dat |> 
  #slice(c(1:3,1344:1346)) |> # test on a few rows
  mutate(allele = map(aa, \(x) paste0(x, sep="",collapse = ""))) |> 
  group_by(gene) |> 
  mutate(aa=(factor(unlist(allele))), allele=as.numeric(factor(unlist(allele), ordered = TRUE))) 

# save file
#save(allele_dat, file="../CAP_data/data/allele_dat.RData")
#load("../CAP_data/data/allele_dat.RData")


## 6. prepare for source attribution methods ##
# which genes only have a single allele?
singles <- allele_dat |> summarise(n=length(unique(allele))) |> filter(n<2) |> pull(gene)

aa_allele_dat <- allele_dat |> 
  ungroup() |> 
  select(LabID, Source, gene, allele) |> 
  pivot_wider(names_from = gene, values_from = allele) |> 
  # remove genes with no variation
  select(!contains(singles))
#save(aa_allele_dat, file="../CAP_data/data/aa_allele_dat.RData")
load("../CAP_data/data/aa_allele_dat.RData")



## 5. calculate hamming distances between alleles ##
aa_distance_matrices <- allele_dat |>
  ungroup() |> 
  # convert allele_dat to list form
  group_split(gene) |> 
  (\(x){set_names(x, nm=map(x, (\(y){paste0(first(y$gene))})))})() |> 
  map(dfun)
genes <- aa_allele_dat |> select(starts_with("CAMP")) |> colnames()
aa_distance_matrices <- aa_distance_matrices[genes]

# save file
#save(aa_distance_matrices, file="../CAP_data/data/aa_distance_matrices.RData")
load("../CAP_data/data/aa_distance_matrices.RData")


