### Simulation study CAP ###

# load libraries
library(tidyverse)
library(boot)

# start with 1 gene, 15 alleles (10 observed, 5 unobserved), 10 individuals, 2 sources
# thinking is that for observed alleles CAP will do no better than either CA or PCO
# but for unobserved alleles CAP will be better, 
# if there is a pattern in CA and strong PCO pattern in a different direction to source then CAP should be better

# generate PCO1 and PCO2 scores
PCO1 <- runif(15, -5, 5)
PCO2 <- runif(15, -1, 1)

# generate probability scores
# choose beta so when we inverse logit we get scores close to zero and close to one
# start with beta=2
beta = 2
probs <- inv.logit(beta * PCO2)

# generate counts for each allele from 10 individuals (size=10)
# need to do this even for the unobserved alleles so we can assign a source to the individuals later
# counts of source one (s1)
count_s1 <- map(probs, ~rbinom(n=1, size=10, p=.)) |> unlist()
# counts of source two (s2)
count_s2 <- 10 - count_s1

# merge into data frame
df <- data.frame(allele = factor(1:15), PCO1, PCO2, p=probs, count_s1, count_s2, observed = c(rep("Y",times=10), rep("N",times=5))) |> 
  mutate(source = ifelse(count_s1>count_s2, "s1", "s2"))

# calculate distances
d <- df |> select(PCO1, PCO2) |> dist()

# plot
df  |> ggplot(aes(x=PCO1, y=PCO2)) + geom_point(aes(colour=source))

# fill out individual row data ie what the original data would actually look like
dat <- df |> select(allele, starts_with("count"), observed) |> 
  rename(s1 = count_s1, s2 = count_s2) |> 
  pivot_longer(cols = starts_with("s"), names_to = "source", values_to = "n") |> 
  uncount(n) |> mutate(across(everything(), as.factor))


## run ranger
source("methods/libs_fns.R")
source("methods/tree_predictions.R")  # to pull out individual tree predictions

# prepare data - add id and split into training and test sets
dat.train <- dat |> rownames_to_column("id") |> filter(observed == "Y") |> select(-observed) |> droplevels() 
dat.test <- dat |> rownames_to_column("id") |> filter(observed == "N") |> select(-observed) |> droplevels()

# 1. ranger with CA scores
# calculate CA scores

## Pull out individual tree predictions
sim_CA0 <- tree_fn(dat.train, dat.test, method="ca0", axes=1, class="source", id="id")
sim_PCO <- tree_fn(dat.train, dat.test, d=d, method="pco", axes=1, class="source", id="id")
sim_CAP <- tree_fn(dat.train, dat.test, d=d, method="cap", k=1, mp=95, axes=1, class="source", id="id")


# 2a. ranger with PCO scores

# 2b. ranger with dist for PCO method to compare with previous (should be nearly the same, perhaps some differences as axes not completely orthogonal)

# 3. ranger with dist and CA scores for CAP method

## store misclassification rates for observed alleles and unobserved alleles


### now repeat 100 times


### now increase number of genes
# what about interaction between genes
# and direction of PCO relationship
# do we need to change how we uncount(n) so it's not the same alleles in the same blocks
