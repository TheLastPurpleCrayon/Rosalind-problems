source("helper_functions.R")

# Introduction to Alternative Splicing
fixed.subsets <- function(n, m) return(sum(chooseZ(n, m:n)) %% 1e6)
#fixed.subsets(n, m)

# HELPER FUNCTION: recursively determine the number of (not necessarily perfect)
# noncrossing matchings given an RNA string
motzkin.rna <- function(start, end) {
  if (end - start + 1 < 2) return(1)
  if (!is.na(memo[start, end])) return(memo[start, end])
  
  minus1 <- motzkin.rna(start+1, end) %% 1e6
  
  first <- vec[start]
  ms <- which(vec[(start+1):end] == bp[first]) + start
  if (length(ms) == 0) return(minus1)
  
  tot <- 0
  for (m in ms) {
    tot <- tot + ((motzkin.rna(start+1, m-1) * motzkin.rna(m+1, end)) %% 1e6)
  }
  result <- (minus1 + tot) %% 1e6
  memo[start, end] <<- result
  return(result)
}

# Motzkin numbers and RNA secondary structures
all.noncrossing.matchings <- function(fasta) {
  parse.fasta.v2(fasta)
  vec <<- strsplit(seqs[1], "")[[1]]
  bp <<- c("A" = "U", "U" = "A", "C" = "G", "G" = "C")
  memo <<- matrix(NA, nchar(seqs[1]), nchar(seqs[1]))
  motzkin.rna(1, nchar(seqs[1])) %% 1e6
}
#all.noncrossing.matchings(fasta)
#tictoc(all.noncrossing.matchings(fasta))