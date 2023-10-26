# Recipe for CA unbiased method

# constructor function for our recipe step
step_ca_unbiased_new <- 
  function(terms, role, trained, ref_scores, options, skip, id) {
    step(
      subclass = "ca_unbiased", 
      terms = terms,
      role = role,
      trained = trained,
      ref_scores = ref_scores,
      options = options,
      skip = skip,
      id = id
    )
  }

# user facing function for our recipe
step_ca_unbiased <- function(
    recipe, 
    ..., 
    role = NA, 
    trained = FALSE, 
    ref_scores = NULL, # scores of levels from the training data
    skip = FALSE,
    options = list(),  # TODO: Add default options here (e.g. for PCO)
    id = rand_id("ca_unbiased")
) {
  
  add_step(
    recipe, 
    step_ca_unbiased_new(
      terms = enquos(...),
      trained = trained,
      role = role,
      ref_scores = ref_scores,
      options = options,
      skip = skip,
      id = id
    )
  )
}

# Prep step

# CA rank (i.e. ranger) scoring function reimplemented
encode_ca_unbiased <- function(x, outcome) {
  #  cat("calling encode ca unbiased... with len(x)=", length(x), "len(outcome)=", length(outcome), "\n")
  #  print(x)
  #  print(outcome)
  x <- droplevels(x)
  if (nlevels(x) < 2) {
    # if we only have one level so we can't do anything
    out <- data.frame(level = levels(x),
                      rank = 1)
    return(out)
  }
  ct <- table(level=x,outcome=outcome)
  P <- ct/rowSums(ct)
  S <- cov.wt(P, wt = rowSums(ct))$cov
  pc1 <- eigen(S)$vectors[,1]   ## PCA of weighted covariance matrix of class probabilities
  
  # compute scores and rank them. New levels will have score 0 but still need ranking
  scores <- P %*% pc1
  ranks <- as.numeric(rank(c(scores, 0), ties.method = "first"))
  list(score = data.frame(level = rownames(scores), rank = ranks[-length(ranks)]),
       new = ranks[length(ranks)])
}

prep.step_ca_unbiased <- function(x, training, info = NULL, ...) {
  # grab the columns we're going to prep
  col_names <- recipes_eval_select(x$terms, training, info)
  
  # grab the outcome column
  outcome_name <- info |> filter(role == "outcome") |> pull(variable)
  if (length(outcome_name) != 1) {
    rlang::abort("One variable with role 'outcome' is required")
  }
  
  # check the column types are what we want: we want factor variables
  # (or categorical variables?)
  check_type(training[, col_names], types = c("character", "factor"))
  check_type(training[, outcome_name], types = c("character", "factor"))
  
  # TODO: Implement this stuff if needed for options to the step
  ## We'll use the names later so make sure they are available
  #  if (x$options$names == FALSE) {
  #    rlang::abort("`names` should be set to TRUE")
  #  }
  
  # e.g. number of PCO axes etc.
  #  if (!any(names(x$options) == "probs")) {
  #    x$options$probs <- (0:100)/100
  #  } else {
  #    x$options$probs <- sort(unique(x$options$probs))
  #  }
  
  # OK, now do the actual CA step on each column
  # This just computes the ranks: the actual data transformation
  # of current variables is done in 'prep' or 'bake' not here
  ref_scores <- purrr::map(
    training[, col_names], 
    \(x) encode_ca_unbiased(x, outcome = training |> pull(outcome_name))
  )
  
  ## Use the constructor function to return the updated object. 
  ## Note that `trained` is now set to TRUE
  
  step_ca_unbiased_new(
    terms = x$terms, 
    trained = TRUE,
    role = x$role, 
    ref_scores = ref_scores,
    options = x$options,
    skip = x$skip,
    id = x$id
  )
}

# Bake step: take our scores and apply them as needed to the columns
apply_unbiased_to_column <- function(x, encoding) {
  
  # create the scoring matrix from the levels of x, and score
  # things we haven't seen as effectively zero
  x <- droplevels(x) # ignore levels we don't have in these data
  
  new_level_rank <- encoding$new
  
  ranks <- data.frame(level = x) |>
    left_join(encoding$score, by="level") |>
    replace_na(list(rank = new_level_rank)) |>
    pull(rank)
  
  ranks
}

bake.step_ca_unbiased <- function(object, new_data, ...) {
  col_names <- names(object$ref_scores)
  check_new_data(col_names, object, new_data)
  
  # iterate over and update our current columns
  for (col_name in col_names) {
    new_data[[col_name]] <- apply_unbiased_to_column(
      x = new_data[[col_name]],
      encoding = object$ref_scores[[col_name]]
    )
  }
  
  # new_data will be a tibble when passed to this function. It should also
  # be a tibble on the way out.
  new_data
}
