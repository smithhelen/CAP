#### Plots to help compare methods ###
# Load libraries and functions
source("methods/libs_fns.R")
library(gridExtra)
library(ggplot2)

## Load results for plotting (from tree predictions)
load("../CAP_Data/results/MC_all_trees.Rdata")

# may need to relabel sources and methods
### now can plot - for individual methods plots use "plots" and "mp_plots" functions; for comparative plots use "methods_plots"
# 1. Bar plots for unique vs non-unique prediction for individual methods
# etc "CA0"    "CA0_CC" "CAP"    "CAP_CC" "PCO"    "PCO_CC"
MC_all_trees |> filter(method=="CA0") |>
    ggplot() +
    geom_col(aes(x=uses_unique, y=prop, fill=prediction)) +
    geom_text(aes(x=uses_unique, y=prop, group=prediction, label=paste0("n=",n)),
              position=position_stack(vjust=0.5), size=3) +
    geom_text(data = MC_all_trees |> filter(method=="CA0") |> group_by(method, uses_unique) |> filter(Source==prediction),
              aes(x=uses_unique, y=1, label=paste0("N=",N)), vjust=-0.5, size=3) +
    labs(title = ("CA0 prediction results"), subtitle="Source",
         y="Proportion of tree predictions\n", x="\nUnique alleles used in prediction",
         fill="Prediction") + 
    theme(text = element_text(size=11)) + 
    facet_wrap(~Source)
MC_all_trees |> filter(method=="CA0") |>
    ggplot() +
    geom_col(aes(x=prediction, y=prop, fill=uses_unique), position="dodge") +
    geom_text(aes(x=prediction, y=prop, group=uses_unique, label=paste0("n=",n)),
              position=position_dodge(width=0.9), vjust=-0.5, size=3) +
    geom_text(aes(x=prediction, y=prop, group=uses_unique, label=paste0("N=",N)),
              position=position_dodge(width=0.9), vjust=1.5, size=3) + 
    labs(title = ("CA0 prediction results"), y="Proportion of tree predictions\n", x="\nPrediction",
         fill="Unique alleles", subtitle="Source") + 
    theme(text = element_text(size=11)) + 
    facet_wrap(~Source)
  list(P1, P2)


# 2. Side by side plots comparing prediction accuracy of different methods
p <- MC_all_trees
p1 <- p |> mutate(pch=ifelse(Source==prediction,21,1)) |> 
    ggplot(aes(x=prediction, y=prop, group=uses_unique)) +
    geom_point(aes(colour=uses_unique, shape=I(pch), bg=uses_unique), size=3) +
    geom_line(aes(colour=uses_unique, lty=uses_unique)) +
    theme_bw(base_size = 14) + 
    scale_y_continuous("Proportion of predictions",expand=c(0.03,0.03),limits=c(0,1)) +
    theme(axis.text.x = element_text(angle = 90), legend.position = "top") +
    facet_grid(Source~method, 
               labeller=labeller(method=c(`CA0`="CA-unbiased",`CA0_CC`="CA_unbiased_CC",`PCO`="PCO",`PCO_CC`="PCO_CC",`CAP`="CAP",`CAP_CC`="CAP_CC"),
                                 Source=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
    scale_colour_manual(values=c("#d34728","#193d87"),labels = c("No", "Yes")) + 
    scale_fill_manual(values=c("#d34728","#193d87"),labels = c("No", "Yes")) + 
    labs(x="Predicted Source", colour="Absent levels used in prediction", bg="Absent levels used in prediction", 
         linetype="Absent levels used in prediction")
p1

p2 <- p |>
  filter(method %in% c("CA0", "PCO", "CAP")) |> 
  mutate(method = factor(method, levels=c("CA0","PCO","CAP"))) |> 
  mutate(pch=ifelse(Source==prediction,21,1)) |> 
  ggplot(aes(x=prediction, y=prop, group=uses_unique)) +
  geom_point(aes(colour=uses_unique, shape=I(pch), bg=uses_unique), size=3) +
  geom_line(aes(colour=uses_unique, lty=uses_unique)) +
  theme_bw(base_size = 14) + 
  scale_y_continuous("Proportion of predictions",expand=c(0.03,0.03),limits=c(0,1)) +
  theme(axis.text.x = element_text(angle = 90), legend.position = "top") +
  facet_grid(Source~method, 
             labeller=labeller(method=c(`CA0`="CA-unbiased",`PCO`="PCO",`CAP`="CAP"),
                               Source=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
  scale_colour_manual(values=c("#d34728","#193d87"),labels = c("No", "Yes")) + 
  scale_fill_manual(values=c("#d34728","#193d87"),labels = c("No", "Yes")) + 
  labs(x="Predicted Source", colour="Absent levels used in prediction", bg="Absent levels used in prediction", 
       linetype="Absent levels used in prediction")
p2

p3 <- p |>
  filter(method %in% c("CA0_CC", "PCO_CC", "CAP_CC")) |> 
  mutate(method = factor(method, levels=c("CA0_CC","PCO_CC","CAP_CC"))) |> 
  mutate(pch=ifelse(Source==prediction,21,1)) |> 
  ggplot(aes(x=prediction, y=prop, group=uses_unique)) +
  geom_point(aes(colour=uses_unique, shape=I(pch), bg=uses_unique), size=3) +
  geom_line(aes(colour=uses_unique, lty=uses_unique)) +
  theme_bw(base_size = 14) + 
  scale_y_continuous("Proportion of predictions",expand=c(0.03,0.03),limits=c(0,1)) +
  theme(axis.text.x = element_text(angle = 90), legend.position = "top") +
  facet_grid(Source~method, 
             labeller=labeller(method=c(`CA0_CC`="CA_unbiased_CC",`PCO_CC`="PCO_CC",`CAP_CC`="CAP_CC"),
                               Source=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
  scale_colour_manual(values=c("#d34728","#193d87"),labels = c("No", "Yes")) + 
  scale_fill_manual(values=c("#d34728","#193d87"),labels = c("No", "Yes")) + 
  labs(x="Predicted Source", colour="Absent levels used in prediction", bg="Absent levels used in prediction", 
       linetype="Absent levels used in prediction")
p3

p4 <- p |>
  mutate(residualised = factor(case_when(grepl("CC",method)~"Yes", TRUE~"No"))) |> 
  mutate(method = factor(case_when(
    grepl("CA0",method)~"CA0",
    grepl("PCO",method)~"PCO",
    grepl("CAP",method)~"CAP"))) |>
  mutate(method = factor(method, levels=c("CA0","PCO","CAP"))) |> 
  mutate(pch=ifelse(Source==prediction,21,1))
p4 |> 
  ggplot(aes(x=prediction, y=prop)) +
  geom_point(aes(colour=method, shape=I(pch), bg=method), size=3) +
  facet_grid(Source~residualised, 
             labeller=labeller(residualised=c(`Yes`="Residualised on CC",`No`="Not residualised"),
                               Source=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep")))
  
  geom_line(aes(colour=method, lty=uses_unique)) +
  theme_bw(base_size = 14) + 
  scale_y_continuous("Proportion of predictions",expand=c(0.03,0.03),limits=c(0,1)) +
  theme(axis.text.x = element_text(angle = 90), legend.position = "top") +
  facet_grid(Source~residualised, 
             labeller=labeller(residualised=c(`Yes`="Residualised on CC",`No`="Not residualised"),
                               Source=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
  scale_colour_manual(values=c("#d34728","#193d87"),labels = c("No", "Yes")) + 
  scale_fill_manual(values=c("#d34728","#193d87"),labels = c("No", "Yes")) + 
  labs(x="Predicted Source", colour="Absent levels used in prediction", bg="Absent levels used in prediction", 
       linetype="Absent levels used in prediction")
p4

p5 <- p  |>
  mutate(residualised = factor(case_when(grepl("CC",method)~"Yes", TRUE~"No"))) |> 
  mutate(method = factor(case_when(
    grepl("CA0",method)~"CA0",
    grepl("PCO",method)~"PCO",
    grepl("CAP",method)~"CAP")))|>
  mutate(Col = case_when(
    method == "CA0" ~ "Dodger blue3",
    method == "PCO" ~ "Orange",
    method == "CAP" ~ "Green4",
    TRUE ~ "Black")) |> 
  filter(Source==prediction) |>
  ggplot(aes(x=uses_unique, y=prop)) + 
  geom_point(aes(colour=method, shape=factor(residualised)),size=3.2, alpha=0.5) + 
  scale_colour_manual(breaks = p2$method, 
                      values = p2$Col) +
  scale_shape_manual(values = c(19,21)) +
  facet_wrap(~Source, nrow = 3) + 
  labs(title = "Accuracy of tree predictions", 
       x="\nUnique alleles used in tree", y="Proportion\n")
p5

pdf(file = "../CAP_Data/results/Results_plots.pdf")
list(p1,p2)
dev.off()



# 3. Side by side plots comparing prediction accuracy of different values of mp
mp_plots(MC_all_trees)

mp_plots <- function(dat){map_dfr(dat, function(x) {
  mp <- x |> map("mp") |> unlist() |> unique()
  answer_all <- x |> map("answer")  |> bind_rows(.id = "Fold") 
  p <- answer_all |> 
    mutate(mp = rep(mp, times = nrow(answer_all))) |>
    mutate(uses_unique = if_else(uses_unique == 0, "No", "Yes")) |>
    mutate(across(Source, factor, levels = c("Poultry","Sheep", "Beef"))) |>
    mutate(across(prediction, factor, levels = c("Poultry","Sheep", "Beef"))) |>
    group_by(uses_unique,Source, prediction,mp) |>
    summarise(n=n()) |>
    group_by(uses_unique,Source) |>
    mutate(N = sum(n)) |>
    mutate(prop = n/N)}) |> ggplot() + 
    geom_col(aes(y=prop, x=Source, fill=prediction)) + 
    labs(title = paste0("CAP prediction results for different mp values"),
         y="Proportion of tree predictions\n", x="\nSource", fill="Prediction") + 
    theme(text = element_text(size=11)) + 
    facet_grid(uses_unique ~ mp, labeller=label_both)
}






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
list(p1,p2,p3,p4,p5,p6,p7)

########################
pdf(file = "../CAP_Data/results/MC_plots.pdf")
list(p1,p2,p3,p4,p5,p6,p7)
dev.off()







