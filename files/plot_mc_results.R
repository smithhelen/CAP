#### Plots to help compare methods ###
# Load libraries and functions
source("methods/libs_fns.R")
library(gridExtra)
#library(ggplot2)

## Load results for plotting 
load("../CAP_Data/results/cgMLST_mc_results.Rdata")

## prepare data for plotting
plot_dat <- cgMLST_mc_results |> 
  mutate(method_long = factor(method, levels=c("CA0","PCO","CAP",
                                               "PCO_Ph", "CAP_Ph",
                                               "CA0_CC","PCO_CC","CAP_CC", 
                                               "CA0_sp","PCO_sp","CAP_sp", 
                                               "CA01","CA02","PCO1","PCO2","CAP1","CAP2",
                                               "CAP_0","CAP_50","CAP_80","CAP_95","CAP_99",
                                               "CA0_aa", "CAP_aa", "PCO_aa",
                                               "CA0_aa_CC", "PCO_aa_CC", "CAP_aa_CC" ,
                                               "CA0_j", "PCO_j", "CAP_j")), .before=1) |> 
  mutate(
    method = substr(method_long, 1, 3) |> factor(levels=c("CA0","PCO","CAP")),
    species = ifelse(method_long |> str_detect("j"),"j","jc") |> factor(levels=c("j", "jc")),
    axes = ifelse(method_long |> str_detect("1"),1,2) |> factor(levels=c("1","2")),
    residualised = ifelse(method_long |> str_detect("CC"),"CC", ifelse(method_long |> str_detect("sp"), "species", "no")) |> factor(levels=c("no", "species", "CC")),
    level = ifelse(method_long |> str_detect("aa"),"aa","nt") |> factor(levels=c("nt", "aa")), 
    mp = ifelse(method_long |> str_extract("_[0-9]+") |> is.na(), 95, method_long |> str_extract("[0-9]+") |> as.numeric()) |> factor(levels = c("0","50","80", "95","99")),  
    .before=2) |> 
  mutate(pch = ifelse(Source==prediction, 21,1)) |> 
  pivot_longer(cols=c(tree.absent, tree.no.absent, forest), names_to = "result", values_to = "p") |> 
  mutate(across(where(is.character), factor)) 





#forest results
p1 <- plot_dat |> filter(method_long %in% c("CA0", "PCO", "CAP")) |> 
  filter(result == "forest") |> 
  filter(prediction == Source) |> 
  ggplot(aes(x=prediction, y=p)) + 
  geom_point(aes(colour=Source, shape=I(pch), bg=Source), size=3) +
  theme_bw(base_size = 11) + 
  scale_y_continuous("Proportion of predictions\n",expand=c(0.03,0.03),limits=c(0,1)) +
  #theme(legend.position = "bottom") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(strip.text.y = element_text(size = 8)) +
  facet_wrap(~method, nrow=1) +
  labs(x="\nPredicted Source", colour="True Source", bg="True Source")

p1

#tree results
p2 <- plot_dat |> filter(method_long %in% c("CA0", "PCO", "CAP")) |> 
  filter(result != "forest") |> 
  ggplot(aes(x=prediction, y=p, group=result)) + 
  geom_point(aes(colour=result, shape=I(pch), bg=result), size=3) +
  geom_line(aes(colour=result)) +
  theme_bw(base_size = 11) + 
  scale_y_continuous("Proportion of predictions\n",expand=c(0.03,0.03),limits=c(0,1)) +
  theme(legend.position = "bottom") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.text.y = element_text(size = 8)) +
  facet_grid(Source~method) +
  scale_colour_manual(values=c("#9d1001","#019d10"), labels=c("Yes", "No")) + 
  scale_fill_manual(values=c("#9d1001","#019d10"), labels=c("Yes", "No")) + 
  labs(x="\nPredicted census region", colour="Absent levels used in prediction", bg="Absent levels used in prediction")

p2

png("plot1.png", width=960, height = 1024)
p1
dev.off()

png("plot2.png", width=960, height = 1024)
p2
dev.off()


## plots
# all methods
p1 <- plot_dat |> 
  ggplot(aes(x=prediction, y=prop, group=uses_unique)) +
  geom_point(aes(colour=uses_unique, shape=I(pch), bg=uses_unique), size=3) +
  geom_line(aes(colour=uses_unique, lty=uses_unique)) +
  theme_bw(base_size = 11) + 
  scale_y_continuous("Proportion of predictions\n",expand=c(0.03,0.03),limits=c(0,1)) +
  theme(axis.text.x = element_text(angle = 90, size=7), 
        legend.position = "bottom", 
        strip.text.x = element_text(angle = 90, size = 8)) +
  facet_grid(Source~method_long, 
             labeller=labeller(Source=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
  scale_colour_manual(values=c("#9d1d01","#003ab3"),labels = c("No", "Yes")) + 
  scale_fill_manual(values=c("#9d1d01","#003ab3"),labels = c("No", "Yes")) + 
  labs(x="\nPredicted Source", colour="Absent levels used in prediction", bg="Absent levels used in prediction", 
       linetype="Absent levels used in prediction", subtitle = "How do the methods compare?")
p1

# CA0 vs PCO vs CAP
p2 <- plot_dat |>
  filter(method_long %in% c("CA0", "PCO", "CAP")) |> 
  mutate(pch=ifelse(Source==prediction,21,1)) |> 
  ggplot(aes(x=prediction, y=prop, group=uses_unique)) +
  geom_point(aes(colour=uses_unique, shape=I(pch), bg=uses_unique), size=3) +
  geom_line(aes(colour=uses_unique, lty=uses_unique)) +
  theme_bw(base_size = 11) + 
  scale_y_continuous("Proportion of predictions\n",expand=c(0.03,0.03),limits=c(0,1)) +
  theme(legend.position = "bottom") +
  facet_grid(Source~method, 
             labeller=labeller(method=c(`CA0`="CA-unbiased",`PCO`="PCO",`CAP`="CAP"),
                               Source=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
  scale_colour_manual(values=c("#9d1d01","#003ab3"),labels = c("No", "Yes")) + 
  scale_fill_manual(values=c("#9d1d01","#003ab3"),labels = c("No", "Yes")) + 
  labs(x="\nPredicted Source", colour="Absent levels used in prediction", bg="Absent levels used in prediction", 
       linetype="Absent levels used in prediction", subtitle = "Comparison of the three base methods")
p2

# CA0 vs PCO vs CAP for CC residualised
p3 <- plot_dat |>
  filter(grepl("CC", method_long), level == "nt") |> 
  ggplot(aes(x=prediction, y=prop, group=uses_unique)) +
  geom_point(aes(colour=uses_unique, shape=I(pch), bg=uses_unique), size=3) +
  geom_line(aes(colour=uses_unique, lty=uses_unique)) +
  theme_bw(base_size = 11) + 
  scale_y_continuous("Proportion of predictions\n",expand=c(0.03,0.03),limits=c(0,1)) +
  theme(legend.position = "bottom") +
  facet_grid(Source~method_long, 
             labeller=labeller(method_long=c(`CA0_CC`="CA_unbiased_CC",`PCO_CC`="PCO_CC",`CAP_CC`="CAP_CC"),
                               Source=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
  scale_colour_manual(values=c("#9d1d01","#003ab3"),labels = c("No", "Yes")) + 
  scale_fill_manual(values=c("#9d1d01","#003ab3"),labels = c("No", "Yes")) + 
  labs(x="\nPredicted Source", colour="Absent levels used in prediction", bg="Absent levels used in prediction", 
       linetype="Absent levels used in prediction", subtitle = "Comparison of the three methods following residualisation (CC)")
p3

# CA0 vs PCO vs CAP for species residualised
p4 <- plot_dat |>
  filter(grepl("sp", method_long), level == "nt") |> 
  ggplot(aes(x=prediction, y=prop, group=uses_unique)) +
  geom_point(aes(colour=uses_unique, shape=I(pch), bg=uses_unique), size=3) +
  geom_line(aes(colour=uses_unique, lty=uses_unique)) +
  theme_bw(base_size = 11) + 
  scale_y_continuous("Proportion of predictions\n",expand=c(0.03,0.03),limits=c(0,1)) +
  theme(legend.position = "bottom") +
  facet_grid(Source~method_long, 
             labeller=labeller(method_long=c(`CA0_sp`="CA_unbiased_sp",`PCO_sp`="PCO_sp",`CAP_sp`="CAP_sp"),
                               Source=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
  scale_colour_manual(values=c("#9d1d01","#003ab3"),labels = c("No", "Yes")) + 
  scale_fill_manual(values=c("#9d1d01","#003ab3"),labels = c("No", "Yes")) + 
  labs(x="\nPredicted Source", colour="Absent levels used in prediction", bg="Absent levels used in prediction", 
       linetype="Absent levels used in prediction", subtitle = "Comparison of the three methods following residualisation (species)")
p4

# effect of residualistion (on CC or species)
p5 <- plot_dat |>
  filter(level == "nt", grepl("CC|sp", residualised) | method_long %in% c("CA0","PCO","CAP")) |> 
  ggplot(aes(x=prediction, y=prop, group=interaction(uses_unique, method))) +
  geom_point(aes(colour=method, shape=I(pch), bg=method), size=3) +
  geom_line(aes(colour=method, lty=uses_unique), linewidth=0.5, alpha=0.8) +
  theme_bw(base_size = 11) + 
  facet_grid(Source~residualised, 
             labeller=labeller(residualised=c(`CC`="Residualised on CC",`species`="Residualised on species",`no`="Not residualised"),
                               Source=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
  scale_y_continuous("Proportion of predictions\n",expand=c(0.03,0.03),limits=c(0,1)) +
  theme(legend.position = "bottom") +
  scale_colour_manual(values=c("#9d1d01","#003ab3", "#f5b200"),labels = c(`CA0`="CA0", `PCO`="PCO", `CAP`="CAP")) + 
  scale_fill_manual(values=c("#9d1d01","#003ab3", "#f5b200"),labels = c(`CA0`="CA0", `PCO`="PCO", `CAP`="CAP")) + 
  labs(x="\nPredicted Source", colour="Method", bg="Method", 
       linetype="Absent levels used in prediction", subtitle = "The effect of residualisation")
p5

# effect of residualistion (on CC)
p6 <- plot_dat  |>
  filter(Source==prediction) |>
  filter(level == "nt", grepl("CC|sp", residualised) | method_long %in% c("CA0","PCO","CAP")) |> 
  ggplot(aes(x=uses_unique, y=prop)) + 
  geom_point(aes(colour=residualised, shape=residualised),size=3.2) + 
  theme_bw(base_size = 11) + 
  scale_colour_manual(values=c("#f5b200", "#003ab3", "#9d1d01"),labels = c(`no`="Not residualised", `species`= "Residualised on species", `CC`="Residualised on CC")) +
  scale_shape_manual(values = c(19,21,4), labels=c(`no`="Not residualised", `species`= "Residualised on species", `CC`="Residualised on CC")) +
  scale_y_continuous("Accuracy\n",expand=c(0.03,0.03),limits=c(0,1)) +
  facet_grid(Source~method, 
             labeller=labeller(Source=c(Cattle="Cattle",Chicken="Chicken",Sheep="Sheep"))) +
  theme(legend.position = "bottom") +
  labs(colour="Residualised", shape="Residualised", 
       x="\nUnique alleles used in tree", subtitle = "The effect of residualisation on overall tree prediction accuracy")
p6

# effect of downweighting recombinant regions
p7 <- plot_dat  |>
  filter(method_long %in% c("CAP","CAP_Ph")) |> 
  ggplot(aes(x=prediction, y=prop, group=uses_unique)) +
  geom_point(aes(colour=uses_unique, shape=I(pch), bg=uses_unique), size=3) +
  geom_line(aes(colour=uses_unique, lty=uses_unique)) +
  theme_bw(base_size = 11) + 
  scale_y_continuous("Proportion of predictions\n",expand=c(0.03,0.03),limits=c(0,1)) +
  theme(legend.position = "bottom") +
  facet_grid(Source~method_long, 
             labeller=labeller(method_long=c(`CAP`="Original Hamming distance", `CAP_Ph`="Phandango adjusted"),
                               Source=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
  scale_colour_manual(values=c("#9d1d01","#003ab3"),labels = c("No", "Yes")) + 
  scale_fill_manual(values=c("#9d1d01","#003ab3"),labels = c("No", "Yes")) + 
  labs(x="\nPredicted Source", colour="Absent levels used in prediction", bg="Absent levels used in prediction", 
       linetype="Absent levels used in prediction", subtitle = "What is the effect of recombination?")
p7

# 1 vs 2 axes 
p8 <- plot_dat |> 
  filter(method_long %in% c("CA01","CA02","PCO1","PCO2","CAP1","CAP2")) |> 
  ggplot(aes(x=prediction, y=prop, group=interaction(factor(axes), uses_unique))) +
  geom_point(aes(colour=factor(axes), shape=I(pch), bg=factor(axes)), size=3, alpha=0.8) +
  geom_line(aes(colour=factor(axes), lty=uses_unique), alpha=0.8) +
  facet_grid(Source~method,
             labeller=labeller(method=c(`CA0`="CA-unbiased",`PCO`="PCO",`CAP`="CAP"),
                               Source=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_y_continuous("Proportion of predictions\n",expand=c(0.03,0.03),limits=c(0,1)) +
  scale_color_manual(values = c("#9d1d01","#003ab3"), label=c("1 axis", "2 axes")) +
  scale_fill_manual(values=c("#9d1d01","#003ab3"), label=c("1 axis", "2 axes")) + 
  labs(x="\nPredicted Source", subtitle = "Does adding a second axis help?", linetype="Absent levels used in prediction", 
       colour = "Dimension", bg="Dimension") 
p8

p8b <- plot_dat2 |> 
  filter(method_long %in% c("CA01","CA02","PCO1","PCO2","CAP1","CAP2")) |> 
  ggplot(aes(x=prediction, y=p, group=axes)) +
  geom_point(aes(colour=factor(axes), shape=I(pch), bg=factor(axes)), size=3, alpha=0.8) +
  geom_line(aes(colour=factor(axes)), alpha=0.8) +
  facet_grid(truth~method,
             labeller=labeller(method=c(`CA0`="CA-unbiased",`PCO`="PCO",`CAP`="CAP"),
                               Source=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_y_continuous("Proportion of predictions\n",expand=c(0.03,0.03),limits=c(0,1)) +
  scale_color_manual(values = c("#9d1d01","#003ab3"), label=c("1 axis", "2 axes")) +
  scale_fill_manual(values=c("#9d1d01","#003ab3"), label=c("1 axis", "2 axes")) + 
  labs(x="\nPredicted Source", subtitle = "Does adding a second axis help?", linetype="Absent levels used in prediction", 
       colour = "Dimension", bg="Dimension") 
p8b

# CAP different mp values 
p9 <- plot_dat |> 
  filter(grepl("_[0-9]+", method_long)) |> 
  ggplot(aes(x=prediction, y=prop, group=interaction(factor(mp), uses_unique))) +
  geom_point(aes(colour=factor(mp), shape=I(pch), bg=factor(mp)), size=3, alpha=0.6) +
  geom_line(aes(colour=factor(mp), lty=uses_unique), alpha=0.6) +
  facet_wrap(~Source,
             labeller=labeller(Source=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
  scale_y_continuous("Proportion of predictions\n",expand=c(0.03,0.03),limits=c(0,1)) +
  scale_color_manual(values =c("red","orange","yellow","green4","blue"), label=c("0", "50","85","95","99")) +
  scale_fill_manual(values=c("red","orange","yellow","green4","blue"), label=c("0", "50","85","95","99")) + 
  labs(x="\nPredicted Source", subtitle = "What is the effect of changing mp (aka PropG)?", linetype="Absent levels used in prediction",
       colour = "mp", bg="mp") +
  theme_bw() +
  theme(legend.position = "bottom")
p9




# plot all (chaos) ###########
p10 <- plot_dat2 |> 
  ggplot(aes(group=method_long, x=prediction, y=p)) +
  geom_point(aes(colour=method_long, shape=I(pch), bg=method_long), size=3) +
  geom_line(aes(colour=method_long, lty=method_long)) +
  facet_wrap(~truth,
             labeller=labeller(truth=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
  scale_y_continuous("Proportion of predictions\n",expand=c(0.03,0.03),limits=c(0,1)) +
  labs(x="\nPredicted Source", y="Proportion of predictions\n", colour = "Method", linetype = "Method", bg="Method") +
  theme(plot.subtitle = element_text(hjust=0.5))
p10


# CAP different mp values ###########
p11 <- plot_dat2 |> 
  ungroup() |> 
  mutate(mp=as.factor(mp)) |> 
  filter(method_long != "CAP_Ph") |> 
  filter(species=="jc" & residualised=="no" & level=="nt" & axes==2 & method=="CAP") |>
  ggplot(aes(x=prediction, y=p)) +
  geom_point(aes(colour=mp, shape=I(pch), bg=mp), size=3, alpha=0.8) +
  geom_line(aes(group=mp, colour=mp), alpha=0.8) +
  facet_wrap(~truth,
             labeller=labeller(truth=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
  scale_y_continuous("Proportion of predictions\n",expand=c(0.03,0.03),limits=c(0,1)) +
  scale_color_manual(values =c("red","orange","yellow","green4","blue"), label=c("0", "50","85","95","99")) +
  scale_fill_manual(values=c("red","orange","yellow","green4","blue"), label=c("0", "50","85","95","99")) + 
  labs(x="\nPredicted Source", subtitle = "What is the effect of changing mp (aka PropG)?", colour = "mp", bg="mp") +
  theme_bw() +
  theme(legend.position = "bottom")
p11

# 3. CA vs PCO vs CAP, all 2 axes ###########
p12 <- plot_dat2 |> 
  filter(method_long != "CAP_Ph") |> 
  filter(species=="jc" & residualised=="no" & level=="nt" & mp==95 & axes==2) |>
  mutate(method = factor(method, levels=c("CA0","PCO","CAP"))) |> 
  ggplot(aes(x=prediction, y=p)) +
  geom_point(aes(colour=method, shape=I(pch), bg=method), size=3, alpha=0.8) +
  geom_line(aes(group=method, colour=method), alpha=0.8) +
  facet_wrap(~truth,
             labeller=labeller(truth=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
  scale_y_continuous("Proportion of predictions\n",expand=c(0.03,0.03),limits=c(0,1)) +
  scale_color_manual(values =c("#9d1d01","#003ab3", "#f5b200"), label=c("CA0","PCO","CAP")) +
  scale_fill_manual(values=c("#9d1d01","#003ab3", "#f5b200"), label=c("CA0","PCO","CAP")) + 
  labs(x="\nPredicted Source", subtitle = "Method comparison (all with 2 axes, mp=95)", colour = "Method", bg="Method") +
  theme_bw() +
  theme(legend.position = "bottom")
p12

# 4. jc vs j species ###########
p13 <- plot_dat2 |> 
  filter(method_long != "CAP_Ph") |> 
  filter(axes==2 & residualised=="no" & level=="nt" & mp==95) |> 
  ggplot(aes(x=prediction, y=p)) +
  geom_point(aes(colour=species, shape=I(pch), bg=species), size=3, alpha=0.8) +
  geom_line(aes(group=species, colour=species), alpha=0.8) +
  facet_grid(truth~method,
             labeller=labeller(method=c(`CA0`="CA-unbiased",`PCO`="PCO",`CAP`="CAP"),
                               truth=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
  scale_y_continuous("Proportion of predictions\n",expand=c(0.03,0.03),limits=c(0,1)) +
  scale_color_manual(values = c("#9d1d01","#003ab3"), label=c("jejuni", "jejuni and coli")) +
  scale_fill_manual(values=c("#9d1d01","#003ab3"), label=c("jejuni", "jejuni and coli")) + 
  labs(x="\nPredicted Source", subtitle = "Does only looking at jejuni help?", colour = "Species", bg="Species") +
  theme_bw() +
  theme(legend.position = "bottom", legend.text = element_text(size=11))
p13

# 5. residualisation CC vs species vs none ###########
p14 <- plot_dat2 |> 
  filter(method_long != "CAP_Ph") |> 
  filter(species=="jc" & axes==2 & level=="nt" & mp==95) |> 
  mutate(residualised = residualised |> factor(levels=c("no","species","CC"))) |> 
  ggplot(aes(x=prediction, y=p)) +
  geom_point(aes(colour=residualised, shape=I(pch), bg=residualised), size=3, alpha=0.8) +
  geom_line(aes(group=residualised, colour=residualised), alpha=0.8) +
  facet_grid(truth~method,
             labeller=labeller(method=c(`CA0`="CA-unbiased",`PCO`="PCO",`CAP`="CAP"),
                               truth=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
  scale_color_manual(values =c("#9d1d01","#003ab3", "#f5b200"), label=c("None","Species","Clonal complex")) +
  scale_fill_manual(values=c("#9d1d01","#003ab3", "#f5b200"), label=c("None","Species","Clonal complex")) + 
  scale_y_continuous("Proportion of predictions\n",expand=c(0.03,0.03),limits=c(0,1)) +
  labs(x="\nPredicted Source", y="Proportion of predictions\n", subtitle = "What happens if we residualise the distance matrices?", colour = "Residualised", bg="Residualised") +
  theme_bw() +
  theme(legend.position = "bottom", legend.text = element_text(size=10))
p14

# 6. aa vs nt level analysis ###########
p15 <- plot_dat2 |> 
  filter(method_long != "CAP_Ph") |> 
  filter(axes==2 & residualised=="no" & mp==95 & species=="jc") |> 
  ggplot(aes(x=prediction, y=p)) +
  geom_point(aes(colour=level, shape=I(pch), bg=level), size=3, alpha=0.8) +
  geom_line(aes(group=level, colour=level), alpha=0.8) +
  facet_grid(truth~method,
             labeller=labeller(method=c(`CA0`="CA-unbiased",`PCO`="PCO",`CAP`="CAP"),
                               truth=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
  scale_color_manual(values = c("#9d1d01","#003ab3"), label=c(`aa`="amino acid", `nt`="nucleotide")) +
  scale_fill_manual(values=c("#9d1d01","#003ab3"), label=c(`aa`="amino acid", `nt`="nucleotide")) + 
  scale_y_continuous("Proportion of predictions\n",expand=c(0.03,0.03),limits=c(0,1)) +
  labs(x="\nPredicted Source", y="Proportion of predictions\n", subtitle = "Does changing the level of analysis help?", colour = "Level", bg="Level") +
  theme_bw() +
  theme(legend.position = "bottom", legend.text = element_text(size=10))
p15

# 7. aa with and without CC residualisation ###########
p16 <- plot_dat2 |> 
  filter(method_long != "CAP_Ph") |> 
  filter(axes==2 & level=="aa" & mp==95 & species=="jc") |> 
  mutate(residualised = residualised |> factor(levels=c("no","species","CC"))) |> 
  ggplot(aes(x=prediction, y=p)) +
  geom_point(aes(colour=residualised, shape=I(pch), bg=residualised), size=3, alpha=0.8) +
  geom_line(aes(group=residualised, colour=residualised), alpha=0.8) +
  facet_grid(truth~method,
             labeller=labeller(method=c(`CA0`="CA-unbiased",`PCO`="PCO",`CAP`="CAP"),
                               truth=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
  scale_color_manual(values = c("#9d1d01","#003ab3"), label=c(`no`="No", `CC`="Clonal complex")) +
  scale_fill_manual(values=c("#9d1d01","#003ab3"), label=c(`no`="No", `CC`="Clonal complex")) + 
  scale_y_continuous("Proportion of predictions\n",expand=c(0.03,0.03),limits=c(0,1)) +
  labs(x="\nPredicted Source", y="Proportion of predictions\n", subtitle = "What about amino acid level with residualisation?", 
       colour = "Residualised", bg="Residualised") +
  theme_bw() +
  theme(legend.position = "bottom", legend.text = element_text(size=10))
p16

# 8. with distance matrices scaled to account for recombination events
p17 <- plot_dat2 |> 
  filter(method_long %in% c("CAP2", "CAP_Ph")) |> 
  mutate(method_long = factor(method_long, levels=c("CAP2","CAP_Ph"))) |> 
  ggplot(aes(x=prediction, y=p))+
  geom_point(aes(colour=method_long, shape=I(pch), bg=method_long), size=3, alpha=0.8) +
  geom_line(aes(group=method_long, colour=method_long), alpha=0.8) +
  facet_wrap(~truth,
             labeller=labeller(truth=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
  scale_color_manual(values = c("#9d1d01","#003ab3"), label=c("Hamming","Phandango")) +
  scale_fill_manual(values=c("#9d1d01","#003ab3"), label=c("Hamming","Phandango")) + 
  scale_y_continuous("Proportion of predictions\n",expand=c(0.03,0.03),limits=c(0,1)) +
  labs(x="\nPredicted Source", y="Proportion of predictions\n", subtitle = "Does accounting for recombination help?", colour = "Distance measure", bg="Distance measure") +
  theme_bw() +
  theme(legend.position = "bottom", legend.text = element_text(size=10))
p17


########################################################################
library(pdftools)
pdf(file = "../CAP_Data/results/Results_plots_1.pdf", width=14, height=10, paper="a4r")
p1
dev.off()

pdf(file = "../CAP_Data/results/Results_plots_2.pdf", width=8, height=10)
list(p10,p2,p12,p3,p4,p5,p14,p9,p11,p7,p17,p15,p16,p8,p8b,p13)
dev.off()

pdf_combine(input = c("../CAP_Data/results/Results_plots_1.pdf", "../CAP_Data/results/Results_plots_2.pdf"), 
            "../CAP_Data/results/Results_plots.pdf")

########################################################################
