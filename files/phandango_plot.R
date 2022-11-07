### Pull out numbers from gff file ###
# load packages
library(ape) # to read gff file
library(tidyverse)

# load data
# start is the first NT of the recombinant region, end is the last NT of the recombinant region, attributes has the number of isolates which are affected
# trial set
dat.gff <- read.gff("C:/Users/user/Dropbox/PhD/SACNZ/Gubbins/Core.full.recombination_predictions.gff") |> select(start, end, attributes)
# full set
dat.gff <- read.gff("C:/Users/user/Dropbox/PhD/SACNZ/Gubbins/allCore.full.recombination_predictions.gff") |> select(start, end, attributes)

# split isolates
dat <- dat.gff |> rownames_to_column("id") |> 
  separate(attributes,sep = "taxa=",into=c("Discard","Isolates")) |> 
  separate(Isolates, sep = "snp_count", into = c("Isolates","Discard2")) |> 
  select(!starts_with("Discard")) |> 
  mutate(Count = str_count(Isolates, "SC")) |> filter(Count > 0)

# split NT position
NT_pos <- function(dat,s,e){
  start = dat |> pull({{s}})
  end = dat |> pull({{e}})
  map2(start, end, ~seq.int(.x,.y)) |> 
    set_names(dat$id)}

# how many times is each position in a recombinant region?
position <- dat |> distinct(start, end, .keep_all = TRUE) |> # remove duplicate rows
  NT_pos(s=start,e=end) |> 
  unlist() |> table() |> as.data.frame() |> mutate(NT = as.numeric(paste0(Var1)))

# plot 1
position |> ggplot(aes(x=NT, y=Freq)) + geom_line()

# scale each position by the number of isolates affected. Need to keep duplicates as different rows had different isolates
position2 <- dat |> NT_pos(s=start,e=end) |> 
  rep(dat$Count) |> unlist() |> table() |> as.data.frame() |> mutate(NT = as.numeric(paste0(Var1)))

# plot 2
position2 |> ggplot(aes(x=NT, y=Freq)) + geom_line()

# save for scaling in random forest
save(position, file="../CAP_data/data/position.Rdata")
save(position2, file="../CAP_data/data/position2.Rdata")
