## How does changing beta effect the outcome?
x1 <- runif(15, -5, 5)
x2 <- runif(15, -1, 1)

df <- data.frame(rep(x2, times=4), beta=rep(c(1,2,20,200), each=15)) |> 
  mutate(x2b = x2 * beta) |> 
  mutate(p=inv.logit(x2b))

group_by(beta) |> arrange(x2b, .by_group = TRUE) |> 
  ungroup() |> mutate(id = rep(1:15, times=4)) 
p1 <- df |> ggplot() + geom_point(aes(x=id, y=p)) + facet_wrap(~beta, ncol=1)
p2 <- df |> ggplot() + geom_point(aes(x=id, y=sort(x2b))) + facet_wrap(~beta, ncol=1)

library(patchwork)

p1_fn <- function(df, b){
  df |> filter(beta==b) |> arrange(p) |>  
    ggplot() + geom_point(aes(x=1:15, y=p)) + labs(title = paste0("beta = ",b), x=NULL) + theme_bw()
}

p2_fn <- function(df, b){
  df |> filter(beta==b) |> arrange(x2b) |> 
    ggplot() + geom_point(aes(x=1:15, y=x2b)) + labs(title = paste0("beta = ",b), y="x", x=NULL) + theme_bw()
}

p <- map(c(1,2,20,200), ~p1_fn(df, .x))
names(p)=c("p1","p2","p3","p4")

px <- map(c(1,2,20,200), ~p2_fn(df, .x))
names(px)=c("p1","p2","p3","p4")

px$p1 + px$p2 + px$p3 + px$p4 + p$p1 + p$p2 + p$p3 + p$p4 + plot_layout(ncol=2, byrow = FALSE)

# generate PCO1 and PCO2 scores
# the direction of greatest variation is along PCO1
PCO1 <- runif(15, -5, 5)
PCO2 <- runif(15, -1, 1)

# choose beta so when we inverse logit we get scores close to zero and close to one
# this is the strength of differentiation along PCO2 between the sources
# we are splitting along PCO2 whereas the direction of greatest variation is along PCO1 
# a higher beta will have a clearer demarcation between the sources
# this demarcation is not the same direction as the greatest spread, so won't be picked up by PCO, but should be picked up by CAP
beta = c(2,20,200)

# generate probability scores - the probability that an allele will belong to source one.
# a higher beta causes probs to group at zero and at one (sigmoid); a lower beta is more linear
probs <- map(beta, ~inv.logit(. * PCO2))

# generate counts for each allele from 10 individuals (size=10)
# based on the probabilities generated from the position on the PCO2 axis, how many will fall in each source
# this is balanced - there are 10 of each allele
count_s1 <- map(probs, function(x) map(x, ~rbinom(n=1, size=10, p=.)) |> unlist())

# fill out individual row data ie what the original data would actually look like
df <- map2(count_s1, probs, function(x,y) data.frame(allele = factor(1:15), PCO1, PCO2, p=y, count_s1=x) )

# plot alleles in PCO space; colour shows likelihood of being in each source
p <- map2(df, c("2","20","200"), function(x,y) x |> ggplot(aes(x=PCO1, y=PCO2)) + geom_point(aes(colour=count_s1))+labs(colour=paste0("beta=",y)))
p[[1]]+p[[2]]+p[[3]]+plot_layout(ncol=1)

# get the CAP scores

dat.train <- dat |> filter(observed == "Y") |> slice_sample(n=75)
dat.test <- dat |> anti_join(dat.train)


