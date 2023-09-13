### Phandango ###

## load libraries
library(tidyverse)
library(stringr)
library(stringi)
library(seqinr)

## function to get reverse complement
rev_comp <- function(sequence)
{return(stri_reverse(chartr("acgtACGT", "tgcaTGCA", sequence)))}

## read data
# cgmlst_ref_map - blast output
dat1 <- read.csv("Phandango/Phandango_data/cgmlst_ref_map.csv") |> select(-file) |> filter(!is.na(closest_allele))

# ref_alignment - aligned version of cgmlst_ref_map
# the alignments are reversed
dat2 <- read.csv("Phandango/Phandango_data/ref_alignment.csv") |> select(-c(file, gene)) |> filter(!is.na(closest_allele))

# ref_mapping - contains the letters in the sequence and where they map to in the aligned version
mapping_dat <- read.csv("Phandango/Phandango_data/ref_mapping.csv") |> mutate(gene = str_split(gene, "[/.]+", simplify = TRUE)[,2]) 

# merge aligned sequence with original columns
aligned_dat <- left_join(dat1, dat2, by=c("allele", "sequence", "closest_allele", "identity")) |> 
  mutate(gene = substr(closest_allele, 1, 8)) |> 
  # replace allele with closest allele
  mutate(allele = str_split(closest_allele, "_", simplify = TRUE)[,2]) |> 
  select(gene, allele, identity, qstart, qend, sstart, send, sequence, aligned_sequence)

# qstart and qend are the positions on the phandango strain where the identified genes are starting and stopping
# sstart and send are the positions on the pubMLST known alleles where the phandango genes are starting and stopping

# there are 22 genes which don't have perfect allele matches
g <- aligned_dat |> filter(identity != 1) |> pull(gene)
g

# there are 6 genes which have gaps in their sequences
aligned_dat |> pull(sequence) |> str_detect("-") |> table()
g1 <- aligned_dat |> slice(which(aligned_dat |> pull(sequence) |> str_detect("-"))) |> pull(gene)
g1

# check if the gaps are also in ref_mapping - they aren't
map(g1, function(x) {mapping_dat |> filter(str_detect(gene, x)) |> pull(letter) |> str_flatten() |> str_detect("-")}) |> unlist()
map(g1, function(x) {aligned_dat |> filter(gene == x) |> pull(sequence) |> str_split("", simplify = TRUE) |> str_detect("-") |> which()})

# are the sequences from phandango or pubMLST?
chq <- aligned_dat |> mutate(qn = abs(qstart-qend)+1, sn = abs(sstart-send)+1) |> 
  mutate(seqn=str_length(sequence), alignedn=str_length(aligned_sequence)) |> 
  mutate(rc_sequence = rev_comp(aligned_sequence)) |> 
  mutate(d1 = stringdist(a=sequence, b=aligned_sequence), d2 = stringdist(a=sequence, b=rc_sequence)) |> 
  mutate(d3 = alignedn-seqn)

# there are 11 where qn doesn't equal sn or seqn
chq |> filter(qn != sn | qn != seqn | sn != seqn) |> view()
g2 <- chq |> filter(qn != sn | qn != seqn | sn != seqn) |> pull(gene)
g2
# these are the 6 genes with gaps plus another 5

# there are 6 (not the same 6 as g above) where there are insertions or deletions in BOTH the sequence column and the aligned/rc column, 
# so the seqn-alignedn != stringdist(a=sequence, b=aligned_sequence)
chq  |> rowwise() |> mutate(d4=min(d1,d2)) |> filter(d3 != d4) |> select(-d4) |> view()
g3 <- chq  |> rowwise() |> mutate(d4=min(d1,d2)) |> filter(d3 != d4) |> select(-d4) |> pull(gene)
g3
x <- chq  |> rowwise() |> mutate(d4=min(d1,d2)) |> filter(d3 != d4) |> select(-d4)
x |> filter(gene=="CAMP1514") |> pull(aligned_sequence) |> str_split(pattern = "") |> as.data.frame() |> write.csv(file="a_seq.csv")
x |> filter(gene=="CAMP1514") |> pull(rc_sequence) |> str_split(pattern = "") |> as.data.frame() |> write.csv(file="rc_seq.csv")
x |> filter(gene=="CAMP1514") |> pull(sequence) |> str_split(pattern = "") |> as.data.frame() |> write.csv(file="seq.csv")

# CAMP0066_13 has aligned length shorter than sequence length - at position 11 there is an "A" in the sequence that is not in the aligned_sequence
chq |> filter(gene=="CAMP0066") |> pull(rc_sequence) |> str_split(pattern = "") |> as.data.frame() |> write.csv(file="rc_seq.csv")
chq |> filter(gene=="CAMP0066") |> pull(sequence) |> str_split(pattern = "") |> as.data.frame() |> write.csv(file="seq.csv")

# CAMP1640_92 has sequence length longer than both qn (qstart-qend+1) and sn (sstart-send+1), and also has three gaps
chq |> filter(gene=="CAMP1640") |> pull(rc_sequence) |> str_split(pattern = "") |> as.data.frame() |> write.csv(file="rc_seq.csv")
chq |> filter(gene=="CAMP1640") |> pull(sequence) |> str_split(pattern = "") |> as.data.frame() |> write.csv(file="seq.csv")

# Other than CAMP1640_92, the sequence length is the longer of either sn or qn
# When qn is longer than sn, qn = sequence length (e.g. CAMP0066_13, CAMP1111_36, CAMP1306_29, CAMP1514_28, CAMP1637_26)
# When qn is shorter than sn, sn = sequence length and there is a gap in the sequence (e.g. CAMP0519_762, CAMP0954_26, CAMP1176_112, CAMP1309_161, CAMP1548_39)
# So I still can't tell if the sequence is the phandango sequence with gaps added where it doesn't match the length of the pubMLST sequence (which would make sense except for CAMP1640_92), or if it is the closest pubMLST allele sequence. Are you please able to check the fasta files for a few of these genes?
  
# check all the abnormal genes
check_these_genes <- chq |> filter(gene %in% c(g, g1,g2,g3)) |> rowwise() |> mutate(stringdiff=min(d1,d2)) |> 
  mutate(aligned_rc = case_when(d2<d1 ~ rc_sequence, TRUE ~ aligned_sequence), .keep="unused") |> 
  mutate(lengthdiff=d3, .keep="unused") |> 
  select(gene, allele, identity, sequence, aligned_rc, qn, sn, seqn, alignedn, lengthdiff, stringdiff)
write.csv(check_these_genes, file="check_these_genes.csv")

# CAMP1484
chq |> filter(gene=="CAMP1484") |> pull(rc_sequence) |> str_split(pattern = "") |> as.data.frame() |> write.csv(file="rc_seq.csv")
chq |> filter(gene=="CAMP1484") |> pull(sequence) |> str_split(pattern = "") |> as.data.frame() |> write.csv(file="seq.csv")
chq |> filter(gene=="CAMP1484") |> pull(aligned_sequence) |> str_split(pattern = "") |> as.data.frame() |> write.csv(file="a_seq.csv")

# CAMP1514
chq |> filter(gene=="CAMP1514") |> pull(rc_sequence) |> str_split(pattern = "") |> as.data.frame() |> write.csv(file="rc_seq.csv")
chq |> filter(gene=="CAMP1514") |> pull(sequence) |> str_split(pattern = "") |> as.data.frame() |> write.csv(file="seq.csv")
chq |> filter(gene=="CAMP1514") |> pull(aligned_sequence) |> str_split(pattern = "") |> as.data.frame() |> write.csv(file="a_seq.csv")

# load new data
newdata <- read.csv("Phandango/Phandango_data/blastn_position_mapping.csv")
CAMP1514 <- newdata |> filter(closest_allele == "CAMP1514_28")

# JM mapping
mapping <- read_csv("Phandango/Phandango_data/alignment_to_reference.csv")
nc_003912 <- read.fasta("Phandango/Phandango_data/NC_003912.fasta")$NC_003912

# Check our mapping. The sequence should directly line up with our q_pos (after flipping to complement as needed)
mapping |>
  mutate(q_base = nc_003912[q_pos]) |>
  group_by(gene) |>
  mutate(q_base = if_else(q_revcomp, seqinr::comp(q_base), q_base)) |>
  filter(aligned_base != q_base)

# load score data
load(file="Phandango/Phandango_data/position.Rdata")
load(file="Phandango/Phandango_data/position2.Rdata")

position$NT |> max()
mapping |> filter(!is.na(q_pos)) |> pull(q_pos) |> max()

