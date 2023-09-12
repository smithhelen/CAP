#### Calculate ngram distances between levels of predictor variable ####

# Load libraries
library(stringdist) # for ngram

# function to create distance (ngram) matrices between alleles for each gene 
# 2*#n-grams shared/(# n-grams in x + # n-grams in y)
ngram_fun <- function(gene, dat, seqdat, ngram=3){  #gene is column name
  alleles <- dat |> pull({gene}) |> unique()
  seqs <- seqdat |> pull({gene}) |> unique() |> as.list()
  names(seqs) <- alleles
  # Get all the qgram counts
  seq_qs <- do.call(what=qgrams, args = c((seqs), q=ngram))
  # pair-wise sum(pmin)/sum() stuff.
  d_ngram <- matrix(0, nrow=nrow(seq_qs), ncol=nrow(seq_qs))
  for (i in 1:nrow(seq_qs)) {
    for (j in i:nrow(seq_qs)) {
      d_ngram[i,j] <- 1-2*sum(pmin(seq_qs[i,], seq_qs[j,]))/sum(seq_qs[i,]+seq_qs[j,])
      d_ngram[j,i] <- d_ngram[i,j]
    }
  }
  colnames(d_ngram) <- alleles
  rownames(d_ngram) <- alleles
  return(d_ngram)
}

# #n-grams shared/# n-grams combined
# This is what Cerda defines in his paper but the above method is probably better
ngram_fun.cerda <- function(gene, dat, seqdat, ngram=3){  #gene is column name
  alleles <- dat |> pull({gene}) |> unique() |> droplevels()
  seqs <- seqdat |> pull({gene}) |> unique() |> droplevels() |> as.list()
  names(seqs) <- alleles
  # Get all the qgram counts
  seq_qs <- do.call(what=qgrams, args = c((seqs), q=ngram)) > 0
  # pair-wise sum(pmin)/sum() stuff.
  d_ngram <- matrix(0, nrow=nrow(seq_qs), ncol=nrow(seq_qs))
  for (i in 1:nrow(seq_qs)) {
    for (j in i:nrow(seq_qs)) {
      d_ngram[i,j] <- sum(seq_qs[i,] & seq_qs[j,])/sum(seq_qs[i,]|seq_qs[j,])
      d_ngram[j,i] <- d_ngram[i,j]
    }
  }
  colnames(d_ngram) <- alleles
  rownames(d_ngram) <- alleles
  return(d_ngram)
}

