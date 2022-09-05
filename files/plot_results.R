#### Plots to help compare methods ###
# Load libraries and functions
source("methods/libs_fns.R")
source("methods/methods_plots.R")

## Load results for plotting (from tree predictions and from cross_validation) - these need to have CAP results added
load("../CAP_Data/results/results_cgMLST_jc.Rdata")  #results_cgMLST_jc
load("../CAP_Data/results/MC_results.Rdata")

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
  
# 4. Overall misclassification
MC_plots(MC_results)