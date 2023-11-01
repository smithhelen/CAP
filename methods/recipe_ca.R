# Recipe for CA method

source('methods/helpers.R')

# constructor function for our recipe step
step_ca_new <- 
  function(terms, role, trained, k, objects, options, skip, id) {
    step(
      subclass = "ca", 
      terms = terms,
      role = role,
      trained = trained,
      k = k,
      objects = objects,
      options = options,
      skip = skip,
      id = id
    )
  }

# user facing function for our recipe
step_ca <- function(
    recipe, 
    ..., 
    role = NA, 
    trained = FALSE, 
    k = NULL, # defaults to number of classes - 1
    objects = NULL, # scores of levels from the training data
    skip = FALSE,
    options = list(),  # TODO: Add default options here (e.g. for PCO)
    id = rand_id("ca")
) {
  
  add_step(
    recipe, 
    step_ca_new(
      terms = enquos(...),
      trained = trained,
      role = role,
      k = k,
      objects = objects,
      options = options,
      skip = skip,
      id = id
    )
  )
}

# Prep step

# CA scoring function (i.e. ranger but with scores not ranks) 
encode_ca <- function(var, outcome, k) {
    cat("calling encode ca ... with len(var)=", length(var), "len(outcome)=", length(outcome), "\n")
  #  print(var)
  #  print(outcome)
  #epsilon <- sqrt(.Machine$double.eps) # in helpers
  
  var <- droplevels(var)
  if (nlevels(var) < 2) {
    # if we only have one level so we can't do anything
    return(null)
    }
  ct <- table(level=var, outcome=outcome)
  if(is.null(k)) {
    k <- ncol(ct) - 1
  }
  P <- ct/rowSums(ct)
  S <- cov.wt(P, wt = rowSums(ct))$cov
  
  eigen_S <- eigen_decomp(S, symmetric = TRUE)   ## PCA of weighted covariance matrix of class probabilities
  # Restrict to a maximum of eigenvectors set by "k" (the number of axes) (default is NULL = ncol(ct)-1)
  nlambdas <- min(sum(eigen_S$values > epsilon), k)
  # principal components
  pc <- eigen_S$vectors
  X <- P %*% pc[, seq_len(nlambdas), drop=FALSE] |> as.data.frame() 
  objects <- list(levels = levels(var),
                  X=X)
  objects
  print(objects)
}

prep.step_ca <- function(x, training, info = NULL, ...) {
  # x is the object from the step_ca function, 
  # training is the training set data (tibble),
  # and info is a tibble that has information on the current set of data eg variable name, type, and role
  
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
  # This just computes the scores: the actual data transformation
  # of current variables is done in 'prep' or 'bake' not here
  objects <- purrr::map(
    training[, col_names], 
    \(var) encode_ca(var = var, outcome = training |> pull(outcome_name), k = x$k)
  )
  
  ## Use the constructor function to return the updated object. 
  ## Note that `trained` is now set to TRUE
  
  step_ca_new(
    terms = x$terms, 
    trained = TRUE,
    role = x$role, 
    k = x$k,
    objects = objects,
    options = x$options,
    skip = x$skip,
    id = x$id
  )
}

# Bake step: take our scores and apply them as needed to the columns
apply_ca_to_column <- function(var, encoding) {
  
  new_level_to_ca <- function(new.var_level, X){
    newX <- X |> map(max) |> as.data.frame() |> mutate(across(everything(), \(x) x+1))
    colnames(newX) <- colnames(X)
    new_var_level_score <- data.frame(level = new.var_level, newX)
    new_var_level_score
  }
  
  # create the scoring matrix from the levels of var, and score
  # things we haven't seen as effectively infinite
  var <- droplevels(var) # ignore levels we don't have in these data
  # Now we figure out which levels are new, and which are the usual
  new_levels <- setdiff(levels(var), encoding$levels)
  new_scores <- map(new_levels, ~new_level_to_ca(new.var_level = {.}, X = encoding$X)) |> list_rbind()
  var_level_score <- bind_rows(data.frame(encoding$X) |> rownames_to_column("level"), new_scores)
  new_cols <- data.frame(level = var) |> left_join(var_level_score, by = "level") |> select(-level)
  new_cols
}

bake.step_ca <- function(object, new_data, ...) {
  # object is the updated step function that has been prepped, 
  # new_data is a tibble of data to be processed.
  col_names <- names(object$objects)
  check_new_data(col_names, object, new_data)
  
  # generate some new names
  new_names <- imap(object$objects, \(x, nm) { paste(nm, "ca", seq_len(ncol(x$X)), sep="_") })
  new_tbl <- tibble::new_tibble(x = list(), nrow=nrow(new_data))
  
  # iterate over and generate our new columns
  for (col_name in col_names) {
    i_col <- new_data[[col_name]]
    i_obj <- object$objects[[col_name]]
    if (!is.null(i_obj)) { # only if we need to include this column...
      new_cols <- apply_ca_to_column(var = i_col, encoding = i_obj)
      new_col_names <- new_names[[col_name]]
      colnames(new_cols) <- new_col_names
      new_tbl[new_col_names] <- new_cols
    }
  }
  
  # new_data will be a tibble when passed to this function. It should also
  # be a tibble on the way out.
  # check the new names and produce our final dataset
  new_tbl <- check_name(new_tbl, new_data, object, names(new_tbl))
  new_data <- bind_cols(new_data, new_tbl)
  new_data <- dplyr::select(new_data, -dplyr::all_of(col_names))
  new_data
}
