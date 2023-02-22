### d to PCO to d ###

# load libraries and functions
source("methods/libs_fns.R")
source("methods/tree_predictions.R")  # to pull out individual tree predictions
library(boot) # for inv.logit function

# set seed
set.seed(24)

# generate PCO1 and PCO2 scores
PCO1 <- runif(15, -5, 5)
PCO2 <- runif(15, -1, 1)

# generate probability scores - the probability that an allele will belong to source one.
probs <- inv.logit(2 * PCO2)

# generate counts for each allele from 10 individuals (size=10)
count_s1 <- map(probs, ~rbinom(n=1, size=10, p=.)) |> unlist()
count_s2 <- 10 - count_s1

PCO <- data.frame(PCO1, PCO2)

# calculate distances
d <- PCO |> dist(diag = TRUE) |> as.matrix()

# re-calculate PCO scores from d for comparison
A <- -0.5 * d^2
B <- dbl_center(A)
eigen_B <- eigen_decomp(B, symmetric=TRUE)
lambdas_B <- eigen_B$values[1:2]
Qo <- eigen_B$vectors
Q <- sweep(Qo[, 1:2, drop=FALSE], 2, sqrt(abs(lambdas_B)), "*") |> as.data.frame()
PCO.d <- data.frame(PCO1 = Q$V1, PCO2 = Q$V2)
cor(PCO,PCO.d)
plot(PCO$PCO1, PCO.d$PCO1)
plot(PCO$PCO2, PCO.d$PCO2)

# recalculate d
d2 <- PCO.d |> dist(diag = TRUE) |> as.matrix()

# are they the same?
all.equal(d, d2)

# the distance matrices are the same even when the PCO axes are not.
# There may be some redundancy in the PCO scores, ie they can get the same distances with different values
# There is not a unique solution to the calculation of PCO scores from d

A <- -0.5 * d2^2
B <- dbl_center(A)
eigen_B <- eigen_decomp(B, symmetric=TRUE)
lambdas_B <- eigen_B$values[1:2]
Qo <- eigen_B$vectors
Q <- sweep(Qo[, 1:2, drop=FALSE], 2, sqrt(abs(lambdas_B)), "*") |> as.data.frame()
PCO.d2 <- data.frame(PCO1 = Q$V1, PCO2 = Q$V2)
all.equal(PCO.d, PCO.d2)

ggplot()+
  geom_point(data=PCO, aes(x=PCO1, y=PCO2), col="red", size=2) + 
  geom_point(data=PCO.d, aes(x=PCO1+0.2, y=PCO2+0.2), col="blue", size=2) + theme_bw()

ggplot()+
  geom_text(data=PCO, aes(x=PCO1, y=PCO2, label=rownames(PCO)), col="red", size=3) + 
  geom_text(data=PCO.d, aes(x=PCO1+0.2, y=PCO2+0.2, label=rownames(PCO.d)), col="blue", size=3) + theme_bw()  + coord_equal()

ggplot()+
  geom_point(data=PCO, aes(x=PCO1, y=PCO2), col="red", size=2) + 
  geom_point(data=PCO.d, aes(x=-PCO1+0.2, y=PCO2-0.24), col="blue", size=2) + theme_bw()

ggplot()+
  geom_text(data=PCO, aes(x=PCO1, y=PCO2, label=rownames(PCO)), col="red", size=3) + 
  geom_text(data=PCO.d, aes(x=-PCO1+0.2, y=PCO2-0.24, label=rownames(PCO.d)), col="blue", size=3) + theme_bw()


one <- PCO[c(1,14),]
two <- PCO.d[c(1,14),]
dist(one)
dist(two)
sqrt(((4.4629345+0.9202889)^2)+((-0.4236312+0.02979553)^2))     
sqrt(((4.261476 +1.100426 )^2)+((-0.8609698+0.2411175)^2))     
ggplot() + geom_point(data=one,aes(x=PCO1,y=PCO2), colour="red") + 
  geom_point(data=two,aes(x=PCO1,y=PCO2), colour="blue") + coord_equal()

# what about cmdscale?
PCO2 <- cmdscale(dist(PCO))
cor(PCO,PCO2)
