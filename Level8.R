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

# Introduction to Set Operations
setOps <- function(string) {
  string <- str_remove_all(string, "\r")
  spl <- strsplit(string, "\n")[[1]]
  n <- as.numeric(spl[1])
  a <- as.numeric(strsplit(substr(spl[2], 2, nchar(spl[2])-1), ",")[[1]])
  b <- as.numeric(strsplit(substr(spl[3], 2, nchar(spl[3])-1), ",")[[1]])
  
  out <- paste0("{", paste(union(a, b), collapse = ", "), "}")
  out <- paste0(out, "\n", paste0("{", paste(intersect(a, b), collapse = ", "), "}"))
  out <- paste0(out, "\n", paste0("{", paste(setdiff(a, b), collapse = ", "), "}"))
  out <- paste0(out, "\n", paste0("{", paste(setdiff(b, a), collapse = ", "), "}"))
  out <- paste0(out, "\n", paste0("{", paste(setdiff(1:n, a), collapse = ", "), "}"))
  out <- paste0(out, "\n", paste0("{", paste(setdiff(1:n, b), collapse = ", "), "}"))
  cat(out)
}
#setOps(string)

# inferring protein from spectrum
spec_to_prot <- function(string){
  chart <- read_delim("useful charts/protein_masses_(monoisotopic).txt", 
                      delim = "   ", col_names = F, show_col_types = F)
  names(chart) <- c("prot", "mass")
  
  spect <- as.numeric(strsplit(string, "\n")[[1]])
  spect <- diff(spect)
  
  out <- character(0)
  for (i in 1:length(spect)) {
    out <- c(out, chart$prot[which.min(abs(spect[i] - chart$mass))])
  }
  cat(out, sep = "")
}
#spec_to_prot(string)
