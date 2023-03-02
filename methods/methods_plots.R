# Plot functions for showing and comparing accuracy of methods

#### Plots to help compare methods ###
# Load libraries and functions
#source("methods/libs_fns.R")
library(gridExtra)
library(ggplot2)

# 1. Bar plots for unique vs non-unique prediction for individual methods
plots <- function(p, method="method") {
  P1 <- p |> filter(method=={{method}}) |>
    ggplot() +
    geom_col(aes(x=uses_unique, y=prop, fill=prediction)) +
    geom_text(aes(x=uses_unique, y=prop, group=prediction, label=paste0("n=",n)),
              position=position_stack(vjust=0.5), size=3) +
    geom_text(data = p |> filter(method=={{method}}) |> group_by(method, uses_unique) |> filter(Source==prediction),
              aes(x=uses_unique, y=1, label=paste0("N=",N)), vjust=-0.5, size=3) +
    labs(title = paste0(method, " prediction results"), subtitle="Source",
         y="Proportion of tree predictions\n", x="\nUnique alleles used in prediction",
         fill="Prediction") + 
    theme(text = element_text(size=11)) + 
    facet_wrap(~Source)
  P2 <- p |> filter(method=={{method}}) |>
    ggplot() +
    geom_col(aes(x=prediction, y=prop, fill=uses_unique), position="dodge") +
    geom_text(aes(x=prediction, y=prop, group=uses_unique, label=paste0("n=",n)),
              position=position_dodge(width=0.9), vjust=-0.5, size=3) +
    geom_text(aes(x=prediction, y=prop, group=uses_unique, label=paste0("N=",N)),
              position=position_dodge(width=0.9), vjust=1.5, size=3) + 
    labs(title = paste0(method, " prediction results"), y="Proportion of tree predictions\n", x="\nPrediction",
         fill="Unique alleles", subtitle="Source") + 
    theme(text = element_text(size=11)) + 
    facet_wrap(~Source)
  list(P1, P2)
}

# 2. Side by side plots comparing prediction accuracy of different methods
methods_plots <- function(p) {
  p1 <- p |> 
    ggplot() + 
    geom_col(aes(x=Source, y=prop, fill=prediction)) + 
    geom_text(aes(x=Source, y=prop, group=prediction, label=paste0("n=",n)),
              position=position_stack(vjust=0.5), size=3) +
    geom_text(data = p |> filter(Source==prediction),
              aes(x=Source, y=1, label=paste0("N=",N)),
              vjust=-0.5, size=3) + 
    labs(title = paste0("Prediction results for different methods"),
         y="Proportion of tree predictions\n", x="\nSource", fill="Prediction") + 
    theme(text = element_text(size=11)) + 
    facet_grid(uses_unique ~ method, labeller=label_both)
  
  p2 <- p |> 
    ggplot() + 
    geom_col(aes(x=Source, y=prop, fill=prediction), colour="black") + 
    labs(title=NULL,
         y="Proportion of tree predictions\n", x="\nSource", fill="Prediction") + 
    theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14)) + 
    facet_grid(uses_unique ~ method, 
               labeller=labeller(method=label_value, uses_unique=c(`No`="No absent levels used in prediction",`Yes`="At least one absent level used in prediction")))
  
  p3 <- p |>
    ggplot() +
    geom_col(aes(x=prediction, y=prop, fill=uses_unique), position="dodge", colour="black") +
    labs(title = "Method", 
         y="Proportion of predictions\n", x="\nPredicted Source", fill="Absent level") + 
    theme(plot.title = element_text(hjust=0.5)) + 
    facet_grid(Source ~ method, 
               labeller=labeller(method=label_value, Source=c(Beef="True Beef",Poultry="True Poultry",Sheep="True Sheep")))
  
  p4 <- p |> mutate(pch=ifelse(Source==prediction,21,1)) |> 
    ggplot(aes(x=prediction, y=prop, group=uses_unique)) +
    geom_point(aes(colour=uses_unique, shape=I(pch), bg=uses_unique), size=3) +
    geom_line(aes(colour=uses_unique, lty=uses_unique)) +
    theme_bw(base_size = 14) + 
    scale_y_continuous("Proportion of predictions",expand=c(0.03,0.03),limits=c(0,1)) +
    theme(axis.text.x = element_text(angle = 90), legend.position = "top") +
    facet_grid(Source~method, 
               labeller=labeller(method=c(`CA_unbiased`="CA-unbiased",`PCO`="PCO",`CAP`="CAP"),
                                 Source=c(Cattle="True Cattle",Chicken="True Chicken",Sheep="True Sheep"))) +
    scale_colour_manual(values=c("#d34728","#193d87"),labels = c("No", "Yes")) + 
    scale_fill_manual(values=c("#d34728","#193d87"),labels = c("No", "Yes")) + 
    labs(x="Predicted Source", colour="Absent levels used in prediction", bg="Absent levels used in prediction", 
         linetype="Absent levels used in prediction")
  
  p5 <- p |> 
    mutate(True=ifelse(Source==prediction,"Y","N")) |> 
    arrange(method, uses_unique, Source, True,prop) |>
    mutate(order = c(1:3)) |> select(c(-n,-N)) |> 
    mutate(order = ifelse(order==1,"A","B")) |>
    group_by(method, uses_unique, Source, order) |> 
    mutate(n=sum(prop)) |> 
    mutate(plot = ifelse(order=="A",1,ifelse(True=="Y",prop,n))) |>
    ungroup() |> select(-c(n, order, True,prop)) |> ggplot() +
    geom_col(aes(x=Source, y=plot, group=method, fill=prediction), position="dodge",colour="black") +
    facet_grid(uses_unique~method, labeller=labeller(uses_unique=c(`No`="No absent levels used in prediction",`Yes`="At least one absent level used in prediction"))) + 
    theme(legend.position = "bottom", plot.subtitle = element_text(hjust=0.5)) + 
    labs(x="\nTrue Source", fill="Predicted Source", y="Proportion of predictions\n") 

  p6 <- p |>
    mutate(Col = case_when(
      method == "CA.zero1" ~ "Yellow",
      method == "CA.zero2" ~ "Purple",
      method == "PCO1" ~ "Orange",
      method == "PCO2" ~ "Dodger blue3",
      method == "CAP2" ~ "Green4",
      method == "CAP2" ~ "Red",
      TRUE ~ "Black"))|>
    filter(Source==prediction) |> 
    mutate(method = factor(method)) |>
    ggplot(aes(x=uses_unique, y=prop, colour=method)) + 
    geom_point(size=2) + 
    scale_colour_manual(breaks = Plotdat$method, 
                        values = as.character(Plotdat$Col)) +
    facet_wrap(~Source, nrow = 3) + 
    labs(title = "Accuracy of tree predictions", 
         x="\nUnique alleles used in tree", y="Proportion\n")
  
  list(p1, p2, p3, p4, p5, p6)
}


# 3. Side by side plots comparing prediction accuracy of different values of mp
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


# 5. attribution results
attribution_plot <- function(dat){ 
  dat |> 
    pivot_longer(-Prediction, names_to = "Method", values_to = "Freq") |> 
    group_by(Method) |> mutate(Proportion = Freq/sum(Freq)) |>
    # reorder so plot looks better
    mutate(Method = factor(Method, levels=c("CA.zero","CA.zero2","PCO","PCO2","CAP","CAP2"))) |>
    ggplot() + geom_col(aes(x=Method, y=Proportion, fill=Prediction)) + 
    labs(title = "Attribution results from different methods using cgMLST data",y="Proportion of isolates\n", x="\nMethod")
}


