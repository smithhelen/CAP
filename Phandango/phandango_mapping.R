## JM mapping

# load libraries
library(tidyverse)
library(seqinr)

mapping <- read_csv("Phandango/Phandango_data/alignment_to_reference.csv")
nc_003912 <- read.fasta("Phandango/Phandango_data/NC_003912.fasta")$NC_003912

# Check our mapping. The sequence should directly line up with our q_pos (after flipping to complement as needed)
mapping |>
  mutate(q_base = nc_003912[q_pos]) |>
  group_by(gene) |>
  mutate(q_base = if_else(q_revcomp, seqinr::comp(q_base), q_base)) |>
  filter(aligned_base != q_base)

# load score data (from plots)
load(file="Phandango/Phandango_data/position.Rdata") # this is the number of times the NT position was in a recombinant event
load(file="Phandango/Phandango_data/position2.Rdata") # this is the number of times as above but multiplied by the number of isolates which had the recombinant event
# I don't know which one we want :-)

position$NT |> max()
mapping |> filter(!is.na(q_pos)) |> pull(q_pos) |> max()

# join frequencies with mapping data and inverse for score
mapping_score <- mapping |> 
  left_join(position |> rename(q_pos = NT, freq = Freq)) |> 
  left_join(position2 |> rename(q_pos = NT, freq.s = Freq)) |> 
  mutate(score = 1/freq, score.s = 1/freq.s)

# what about scaling freq.s to the same range as freq

# save
save(mapping_score, file="Phandango/Phandango_data/mapping_score.Rdata")
load(file="Phandango/Phandango_data/mapping_score.Rdata") # (mapping_score)

