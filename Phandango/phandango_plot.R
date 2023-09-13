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
# gives a list of every nt position between start and end which can then be tallied
NT_pos <- function(dat,s,e){
  start = dat |> pull({{s}})
  end = dat |> pull({{e}})
  map2(start, end, ~seq.int(.x,.y)) |> 
    set_names(dat$id)}

# how many times is each position in a recombinant region?
position <- dat |> distinct(start, end, .keep_all = TRUE) |> # remove duplicate rows - some isolates are on separate rows for some reason
  NT_pos(s=start,e=end) |> 
  unlist() |> table() |> as.data.frame() |> mutate(NT = as.numeric(paste0(Var1)), .keep="unused", .before = 1)

# plot 1
p1 <- position |> ggplot(aes(x=NT, y=Freq)) + geom_line(col="lightblue") + #geom_line(col="blue") + 
  labs(title = "Phandango - number of times an NT position was in a recombinant event", 
       x="q position", y="Frequency")

# scale each position by the number of isolates affected. Need to keep duplicates as different rows had different isolates
position2 <- dat |> NT_pos(s=start,e=end) |> 
  rep(dat$Count) |> unlist() |> table() |> as.data.frame() |> mutate(NT = as.numeric(paste0(Var1)), .keep="unused", .before = 1)

# plot 2
p2 <- position2 |> ggplot(aes(x=NT, y=Freq)) + geom_line(col="pink") + #geom_line(col="darkgreen") + 
  labs(title = "Phandango - number of times an NT position was in a recombinant event multiplied by the number of isolates affected", 
       x="q position", y="Frequency") 

png("Phandango/Phandango_data/PhandangoPlot.png", width=800, height=640)
gridExtra::grid.arrange(ncol=1, p1, p2)
dev.off()

# save for scaling in random forest
#save(position, file="Phandango/Phandango_data/position.Rdata")
#save(position2, file="Phandango/Phandango_data/position2.Rdata")

# load
load(file="Phandango/Phandango_data/position.Rdata")
load(file="Phandango/Phandango_data/position2.Rdata")

# now join with the mapping data