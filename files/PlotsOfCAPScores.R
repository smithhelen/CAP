### if there is a pattern in CA and strong PCO pattern in a different direction to source then CAP should be better

# load libraries and functions
source("files/simulation_study.R")
library(patchwork)
library(Manu)

dat <- gen_data(beta=200)
d <- list(dat$d)
preds <- run_ranger(dat = dat$dat, ntrees=5, d=d, k=2, m=2, axes=2, class="source", id="id", var_id="allele")

# pco scores
preds$dat.scores$pco$train$extra$allele$Q

# Q is the same as Qo * sqrt(lambdas_B)
preds$dat.scores$cap$train$extra$allele$Qo *
sqrt(preds$dat.scores$pco$train$extra$allele$lambdas_B)

# cap scores
preds$dat.scores$cap$train$extra$allele$C_score
#C_score <- Qo %*% U
preds$dat.scores$cap$train$extra$allele$Qo %*%
  preds$dat.scores$cap$train$extra$allele$U

# plot pco
dat$df |> ggplot(aes(x=PCO1, y=PCO2)) + geom_point(aes(colour=count_s1))+theme_bw() + 
  scale_colour_gradientn(colours = get_pal("Pohutukawa")[c(1,4)])

(preds$dat.scores$cap$train$extra$allele$Qo *
    sqrt(preds$dat.scores$pco$train$extra$allele$lambdas_B))|> as.data.frame() |> rename(PCO1=V1, PCO2=V2) |> rownames_to_column("allele") |> 
  left_join(dat$df |> slice(1:10)|> select(allele, count_s1)) |> 
  ggplot(aes(x=PCO1, y=PCO2)) + geom_point(aes(colour=count_s1))+theme_bw()+
  scale_colour_gradientn(colours = get_pal("Pohutukawa")[c(1,4)])

pco_plot <- preds$dat.scores$pco$train$extra$allele$Q |> as.data.frame() |> rename(PCO1=V1) |> rownames_to_column("allele") |> 
  left_join(dat$df |> slice(1:10) |> select(allele, count_s1)) |> 
  ggplot(aes(x=PCO1, y=c(rep(0, times=10)))) + geom_point(aes(colour=count_s1), size=4)+theme_bw()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+ 
  scale_colour_gradientn(colours = get_pal("Pohutukawa")[c(1,4)])

# plot cap
cap_plot <- preds$dat.scores$cap$train$extra$allele$C_score |> as.data.frame() |> rename(CAP1=V1) |> rownames_to_column("allele") |> 
  left_join(dat$df |> slice(1:10)) |> 
  ggplot(aes(x=CAP1, y=c(rep(0, times=10)))) + geom_point(aes(colour=count_s1), size=4)+theme_bw()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())+ 
  scale_colour_gradientn(colours = get_pal("Pohutukawa")[c(1,4)])

p <- pco_plot + cap_plot + plot_layout(ncol=1)
png("PCO_vs_CAP.png", width=640, height=480)
p
dev.off()
