# load libraries and functions
source("files/simulation_study.R")

### intermediary checks
sim_1 <- sim_fn(beta = 2, ntrees = 5)

# plot alleles in PCO space; colour shows likelihood of being in each source
sim_1$mc.extras$dat$df |> ggplot(aes(x=PCO1, y=PCO2)) + geom_point(aes(colour=count_s1))
sim_1$mc.extras$dat$PCO
sim_1$mc.extras$dat$Q
sim_1$mc.extras$tree_preds
sim_1$mc.extras$dat.scores$pco$train$extra$allele$Q
sim_1$mc
sim_100[[1]]$mc.extras$dat$df


### now repeat 100 times
#sim_100_b2 <- map(1:100, ~sim_fn(beta=2, ntrees=500))
#sim_100_b20 <- map(1:100, ~sim_fn(beta=20, ntrees=500))
#sim_100 <- list(beta_2 = sim_100_b2, beta_20 = sim_100_b20)
#save(sim_100, file="sim_100.Rdata")

#mc_100_b2 <- map(sim_100_b2,"mc") |> bind_rows(.id="simulation") |> mutate(beta = 2)
#mc_100_b20 <- map(sim_100_b20,"mc") |> bind_rows(.id="simulation") |> mutate(beta = 20)
#mc_100 <- bind_rows(mc_100_b2, mc_100_b20)

# and plot
library(Manu)
library(wesanderson)

sim_plot_1 <- mc_100 |> mutate(method = factor(method, levels = c("ca0","pco","pco2","cap"))) |> 
  mutate(method = fct_recode(as.factor(method), CA_unbiased="ca0", PCO="pco", PCO_scores="pco2", CAP="cap")) |> 
  mutate(beta = factor(beta)) |> 
  ggplot() + geom_boxplot(aes(x=uses_unique, y=mc, fill=uses_unique)) + 
  facet_grid(beta~method, labeller = labeller(beta = c(`2`="Beta = 2",`20`="Beta = 20"))) + 
  theme_bw(base_size = 11) + 
  labs(y="Misclassification rate") + 
  theme(legend.position = "top", legend.title = element_blank(),
        axis.ticks.x = element_blank(), axis.title.x = element_blank(),axis.text.x = element_blank()) +
  scale_fill_manual(breaks = c("No","Yes"), 
                    values=manu_palettes$Kaka, #values = c("#d34728","#3b4da1"),
                    labels = c("No absent levels used for prediction", "At least one absent level used for prediction"))

sim_plot_2 <- mc_100 |> mutate(method = factor(method, levels = c("ca0","pco","pco2","cap"))) |> 
  mutate(method = fct_recode(as.factor(method), CA_unbiased="ca0", PCO="pco", PCO_scores="pco2", CAP="cap")) |> 
  mutate(uses_unique = factor(uses_unique, levels=c("No","Yes"))) |> 
  mutate(beta = factor(beta)) |> 
  ggplot() + geom_boxplot(aes(x=method, y=mc, fill=method)) + 
  facet_grid(beta~uses_unique, 
             labeller=labeller(uses_unique=c(`No`="No absent levels used in prediction",`Yes`="At least one absent level used in prediction"),
                               beta = c(`2`="Beta = 2",`20`="Beta = 20"))) + 
  theme_bw(base_size = 11) + 
  labs(y="Misclassification rate", fill="Method") + 
  theme(axis.ticks.x = element_blank(), axis.title.x = element_blank(),axis.text.x = element_blank()) +
  scale_fill_manual(values=manu_palettes$Kaka) #values=wes_palette("FantasticFox1"))

# Save plot
mods <- ls(pattern='sim_plot_[0-9]+[b]?')
#save(list = mods, file = "../CAP_data/results/simulation_plots.RData")

# Load plots ready for png output
load(file = "../CAP_data/results/simulation_plots.RData")

png("SimPlot1.png", width=640, height=480)
sim_plot_1
dev.off()
png("SimPlot2.png", width=640, height=480)
sim_plot_2
dev.off()
