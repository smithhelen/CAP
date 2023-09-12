# plots for simulations

library(Manu)
library(wesanderson)

sim_plot_1 <- function(plotdat){
  plotdat |> mutate(method = factor(method, levels = c("ca0","pco","pco2","cap"))) |> 
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
}

sim_plot_2 <- function(plotdat){
  plotdat |> mutate(method = factor(method, levels = c("ca0","pco","pco2","cap"))) |> 
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
    scale_fill_manual(values=manu_palettes$Kaka) 
}
