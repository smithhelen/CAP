# load libraries and functions
# single gene:
source("files/simulation_study.R")
source("files/plots_for_simulations.R")

### intermediary checks
sim_1 <- sim_fn(beta = 2, ntrees = 5)
sim_10 <- map(1:10, ~sim_fn(beta=2, ntrees=5))
mc_10 <- map(sim_10,"mc") |> bind_rows(.id="simulation") |> mutate(beta = 2)

### now repeat 100 times
#sim_100_b2 <- map(1:100, ~sim_fn(beta=2, ntrees=500))
#sim_100_b20 <- map(1:100, ~sim_fn(beta=20, ntrees=500))
#sim_100 <- list(beta_2 = sim_100_b2, beta_20 = sim_100_b20)
#save(sim_100, file="sim_100.Rdata")

#mc_100_b2 <- map(sim_100_b2,"mc") |> bind_rows(.id="simulation") |> mutate(beta = 2)
#mc_100_b20 <- map(sim_100_b20,"mc") |> bind_rows(.id="simulation") |> mutate(beta = 20)
#mc_100 <- bind_rows(mc_100_b2, mc_100_b20)

# and plot
p1 <- sim_plot_1(mc_100)
p2 <- sim_plot_2(mc_100)

# Save plot
mods <- ls(pattern='sim_plot_[0-9]+[b]?')
#save(list = mods, file = "../CAP_data/results/simulation_plots.RData")

# Load plots ready for png output
load(file = "../CAP_data/results/simulation_plots.RData")

png("SimPlot1.png", width=640, height=480)
p1
dev.off()
png("SimPlot2.png", width=640, height=480)
p2
dev.off()

#######################################################################################################
# two genes:
source("files/simulation_study_2genes.R")

### intermediary checks
sim_1_2genes <- sim_fn_2genes(beta = 2, ntrees = 5)
sim_10_2genes <- map(1:10, ~sim_fn_2genes(beta=2, ntrees=5))
mc_10_2genes <- map(sim_10_2genes,"mc") |> bind_rows(.id="simulation") |> mutate(beta = 2)

sim_plot_1(mc_10_2genes)

### now repeat 100 times
#sim_100_b2_2genes <- map(1:100, ~sim_fn_2genes(beta=2, ntrees=500))
#sim_100_b20_2genes <- map(1:100, ~sim_fn_2genes(beta=20, ntrees=500))
#sim_100_2genes <- list(beta_2 = sim_100_b2, beta_20 = sim_100_b20)
#save(sim_100_2genes, file="sim_100_2genes.Rdata")

#mc_100_b2_2genes <- map(sim_100_b2_2genes,"mc") |> bind_rows(.id="simulation") |> mutate(beta = 2)
#mc_100_b20_2genes <- map(sim_100_b20_2genes,"mc") |> bind_rows(.id="simulation") |> mutate(beta = 20)
#mc_100_2genes <- bind_rows(mc_100_b2_2genes, mc_100_b20_2genes)

# and plot
p1 <- sim_plot_1(mc_100_b2_2genes )
p2 <- sim_plot_2(mc_100_b2_2genes )

# Save plot
#mods <- ls(pattern='sim_plot_[0-9]+[b]?')
#save(list = mods, file = "../CAP_data/results/simulation_plots.RData")

# Load plots ready for png output
#load(file = "../CAP_data/results/simulation_plots.RData")

png("SimPlot1.png", width=640, height=480)
p1
dev.off()
png("SimPlot2.png", width=640, height=480)
p2
dev.off()

