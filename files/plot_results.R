#### Plots to help compare methods ###
# Load libraries and functions
source("methods/libs_fns.R")
source("methods/methods_plots.R")

## Load results for plotting (from tree predictions and from cross_validation) - these need to have CAP results added
load("../CAP_Data/results/results_cgMLST_jc.Rdata")  #results_cgMLST_jc
load("../CAP_Data/results/all_MC_results.Rdata")

# to check differences between CAP with different mp values, will need to recreate results_cgMLST_jc (and for MLST data) 
# in tree_data.R as the above are with mp=95

# calculate tree misclassifications to plot tree data
MC_all_trees <- misclass_tree_fn(results_cgMLST_jc) 

# may need to relabel sources and methods
### now can plot - for individual methods plots use "plots" and "mp_plots" functions; for comparative plots use "methods_plots"
# 1. Bar plots for unique vs non-unique prediction for individual methods
plots(MC_all_trees, method="CA.zero1")
plots(MC_all_trees, method="CA.zero2")
plots(MC_all_trees, method="PCO1")
# etc

# 2. Side by side plots comparing prediction accuracy of different methods - will return list with 6 plots!
methods_plots(MC_all_trees)

# 3. Side by side plots comparing prediction accuracy of different values of mp
mp_plots(MC_all_trees)
  
# 4. Overall misclassification (split by source) # run 3/11/2022
# load results from cross_validation.R
load("../CAP_Data/results/all_MC_results.Rdata")  # from cross_validation.R

# prepare data for plotting
plot_dat <- all_MC_results |> map_dfr("conf.av") |> 
  mutate(method_long = names(all_MC_results),.before=1) |> as.data.frame() |> 
  filter(!method_long %in% c("CAP_95","CAP_CC_95")) |> 
  mutate(
    method = substr(method_long, 1, 3) |> factor(levels=c("CA0","PCO","CAP")),
    species = ifelse(method_long |> str_detect("j"),"j","jc"),
    axes = ifelse(method_long |> str_detect("1"),1,2),
    residualised = ifelse(method_long |> str_detect("CC"),"CC",ifelse(method_long |> str_detect("sp"), "species","no")),
    level = ifelse(method_long |> str_detect("aa"),"aa","nt"), 
    mp = ifelse(method_long |> str_extract("_[0-9]+") |> is.na(), 95, method_long |> str_extract("[0-9]+") |> as.numeric()),  
    .before=2) |> 
  pivot_longer(starts_with(c("Cattle","Chicken","Sheep")), names_to = "T_P", values_to = "p") |> rowwise() |> 
  mutate(truth = str_split(T_P, "_|\\.")[[1]][1], prediction = str_split(T_P, "_|\\.")[[1]][2]) |> 
  mutate(pch = ifelse(truth==prediction, 21,1))

# plot all (chaos) ###########
plot_dat |> 
  ggplot(aes(group=method_long, x=prediction, y=p)) +
  geom_point(aes(colour=method_long, shape=I(pch), bg=method_long), size=3) +
  geom_line(aes(colour=method_long, lty=method_long)) +
  facet_wrap(~truth,
             labeller=labeller(truth=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
  labs(x="\nPredicted Source", y="Proportion of predictions\n", colour = "Method", linetype = "Method", bg="Method") +
    theme(plot.subtitle = element_text(hjust=0.5))
  
# 1. 1 vs 2 axes ###########
p1 <- plot_dat |> 
  filter(species=="jc" & residualised=="no" & level=="nt" & mp==95) |> 
  ggplot(aes(x=prediction, y=p, alpha=0.1)) +
  scale_color_manual(values = c("#d34728","#193d87"), label=c("1 axis", "2 axes")) +
  scale_fill_manual(values=c("#d34728","#193d87"), label=c("1 axis", "2 axes")) + 
  geom_point(aes(colour=factor(axes), shape=I(pch), bg=factor(axes)), size=3, alpha=0.6) +
  geom_line(aes(group=factor(axes), colour=factor(axes)), alpha=0.6) +
  facet_grid(method~truth,
             labeller=labeller(method=c(`CA0`="CA-unbiased",`PCO`="PCO",`CAP`="CAP"),
                               truth=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
   labs(x="\nPredicted Source", y="Proportion of predictions\n", subtitle = "Does adding a second axis help?", colour = "Dimension", bg="Dimension") +
  theme(plot.subtitle = element_text(hjust=0.5)) +
  theme_bw()

# 2. CAP different mp values ###########
p2 <- plot_dat |> ungroup() |> mutate(mp=as.factor(mp)) |> 
  filter(species=="jc" & residualised=="no" & level=="nt" & axes==2 & method=="CAP") |>
  ggplot(aes(x=prediction, y=p, alpha=0.1)) +
  geom_point(aes(colour=mp, shape=I(pch), bg=mp), size=3, alpha=0.6) +
  geom_line(aes(group=mp, colour=mp), alpha=0.6) +
  scale_color_manual(values =c("red","orange","yellow","green4","blue"), label=c("0", "50","85","95","99")) +
  scale_fill_manual(values=c("red","orange","yellow","green4","blue"), label=c("0", "50","85","95","99")) + 
  facet_wrap(~truth,
             labeller=labeller(truth=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
  labs(x="\nPredicted Source", y="Proportion of predictions\n", subtitle = "What is the effect of changing mp (aka PropG)?", colour = "mp", bg="mp") +
  theme(plot.subtitle = element_text(hjust=0.5)) +
  theme_bw()

# 3. CA vs PCO vs CAP, all 2 axes ###########
p3 <- plot_dat |> filter(species=="jc" & residualised=="no" & level=="nt" & mp==95 & axes==2) |>
  ggplot(aes(x=prediction, y=p, alpha=0.1)) +
  geom_point(aes(colour=method, shape=I(pch), bg=method), size=3, alpha=0.6) +
  geom_line(aes(group=method, colour=method), alpha=0.6) +
  scale_color_manual(values =c("red","#21a334","#0048ba"), label=c("CA0","PCO","CAP")) +
  scale_fill_manual(values=c("red","#21a334","#0048ba"), label=c("CA0","PCO","CAP")) + 
  facet_wrap(~truth,
             labeller=labeller(truth=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
  labs(x="\nPredicted Source", y="Proportion of predictions\n", subtitle = "Method comparison (all with 2 axes, mp=95)", colour = "Method", bg="Method") +
  theme(plot.subtitle = element_text(hjust=0.5)) +
  theme_bw()
  
# 4. jc vs j species ###########
p4 <- plot_dat |> 
  filter(axes==2 & residualised=="no" & level=="nt" & mp==95) |> 
  ggplot(aes(x=prediction, y=p, alpha=0.1)) +
  scale_color_manual(values = c("#d34728","#193d87"), label=c("jejuni", "jejuni and coli")) +
  scale_fill_manual(values=c("#d34728","#193d87"), label=c("jejuni", "jejuni and coli")) + 
  geom_point(aes(colour=species, shape=I(pch), bg=species), size=3, alpha=0.6) +
  geom_line(aes(group=species, colour=species), alpha=0.6) +
  facet_grid(method~truth,
             labeller=labeller(method=c(`CA0`="CA-unbiased",`PCO`="PCO",`CAP`="CAP"),
                               truth=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
  labs(x="\nPredicted Source", y="Proportion of predictions\n", subtitle = "Does only looking at jejuni help?", colour = "Species", bg="Species") +
  theme(plot.subtitle = element_text(hjust=0.5)) +
  theme_bw()

# 5. residualisation CC vs species vs none ###########
p5 <- plot_dat |> 
  filter(species=="jc" & axes==2 & level=="nt" & mp==95) |> 
  mutate(residualised = residualised |> factor(levels=c("no","species","CC"))) |> 
  ggplot(aes(x=prediction, y=p, alpha=0.1)) +
  scale_color_manual(values =c("red","#21a334","#0048ba"), label=c("No","Species","Clonal complex")) +
  scale_fill_manual(values=c("red","#21a334","#0048ba"), label=c("No","Species","Clonal complex")) + 
  geom_point(aes(colour=residualised, shape=I(pch), bg=residualised), size=3, alpha=0.6) +
  geom_line(aes(group=residualised, colour=residualised), alpha=0.6) +
  facet_grid(method~truth,
             labeller=labeller(method=c(`CA0`="CA-unbiased",`PCO`="PCO",`CAP`="CAP"),
                               truth=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
  labs(x="\nPredicted Source", y="Proportion of predictions\n", subtitle = "What happens if we residualise the distance matrices?", colour = "Residualised", bg="Residualised") +
  theme(plot.subtitle = element_text(hjust=0.5)) +
  theme_bw()

# 6. aa vs nt level analysis ###########
p6 <- plot_dat |> 
  filter(axes==2 & residualised=="no" & mp==95 & species=="jc") |> 
  ggplot(aes(x=prediction, y=p, alpha=0.1)) +
  scale_color_manual(values = c("#d34728","#193d87"), label=c("amino acid", "nucleotide")) +
  scale_fill_manual(values=c("#d34728","#193d87"), label=c("amino acid", "nucleotide")) + 
  geom_point(aes(colour=level, shape=I(pch), bg=level), size=3, alpha=0.6) +
  geom_line(aes(group=level, colour=level), alpha=0.6) +
  facet_grid(method~truth,
             labeller=labeller(method=c(`CA0`="CA-unbiased",`PCO`="PCO",`CAP`="CAP"),
                               truth=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
  labs(x="\nPredicted Source", y="Proportion of predictions\n", subtitle = "Does changing the level of analysis help?", colour = "Level", bg="Level") +
  theme(plot.subtitle = element_text(hjust=0.5)) +
  theme_bw()

# 7. aa with and without CC residualisation ###########
p7 <- plot_dat |> 
  filter(axes==2 & level=="aa" & mp==95 & species=="jc") |> 
  mutate(residualised = residualised |> factor(levels=c("no","species","CC"))) |> 
  ggplot(aes(x=prediction, y=p)) +
  scale_color_manual(values = c("#d34728","#193d87"), label=c("No", "Clonal complex")) +
  scale_fill_manual(values=c("#d34728","#193d87"), label=c("No", "Clonal complex")) + 
  geom_point(aes(colour=residualised, shape=I(pch), bg=residualised), size=3, alpha=0.6) +
  geom_line(aes(group=residualised, colour=residualised), alpha=0.6) +
  facet_grid(method~truth,
             labeller=labeller(method=c(`CA0`="CA-unbiased",`PCO`="PCO",`CAP`="CAP"),
                               truth=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
  labs(x="\nPredicted Source", y="Proportion of predictions\n", subtitle = "What about amino acid level with residualisation?", 
       colour = "Residualised", bg="Residualised") +
  theme(plot.subtitle = element_text(hjust=0.5)) +
  theme_bw()

########################
pdf(file = "../CAP_Data/results/MC_plots.pdf")
list(p1,p2,p3,p4,p5,p6,p7)
dev.off()







