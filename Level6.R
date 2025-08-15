source("helper_functions.R")

# finding a spliced motif
subseq.idx <- function(fasta) {
  spl <- str_split(fasta, pattern = ">")[[1]] |>
    str_remove_all("\\n") |> 
    str_extract("(?<=Rosalind_\\w{4})\\w*")
  dna <- spl[2]
  subsequence <- str_split(spl[3], "")[[1]]
  
  vec <- numeric(0)
  for (i in 1:length(subsequence)) {
    nextmatch <- str_locate(dna, (subsequence[i]))[1,1] |> as.numeric()
    offset <- if (is.infinite(max(vec))) 0 else max(vec)
    vec <- c(vec, nextmatch+offset)
    dna <- str_sub(dna, start = nextmatch+1, end = str_length(dna))
  }
  str_flatten(vec, collapse = " ")
}
#subseq.idx(fasta)

# transitions and transversions
transition.transversion.ratio <- function(fasta) {
  parse.fasta(fasta)
  
  gene1 <- strsplit(seqs[1], "")[[1]]
  gene2 <- strsplit(seqs[2], "")[[1]]
  
  mut1 <- mut2 <- character(0)
  nTransitions <- nTransversions <- 0
  pur <- c("A", "G") # A and G are purines
  pyr <- c("C", "T") # C ant T are pyrimidines
  for (i in 1:length(gene1)) {
    if (gene1[i] != gene2[i]) {
      mut1 <- c(mut1, gene1[i])
      mut2 <- c(mut2, gene2[i])
    }
  }
  for (j in 1:length(mut1)) {
    if (mut1[j] %in% pur & mut2[j] %in% pur) {
      nTransitions = nTransitions + 1
    } else if (mut1[j] %in% pyr & mut2[j] %in% pyr) {
      nTransitions = nTransitions + 1
    } else {
      nTransversions = nTransversions + 1
    }
  }
  nTransitions/nTransversions
}
#transition.transversion.ratio(fasta)

# enumerating k-mers lexicographically
enumerkate <- function(alphabet, n) {
  alph <- str_remove_all(alphabet, "\\s")
  
  list <- list(strsplit(alph, "")[[1]])
  for (i in 1:(n-1)){
    for (j in length(list[[1]]):1) {
      list[[j]] <- paste0(list[[1]][j], strsplit(alph, "")[[1]])
    }
    list <- list(unlist(list))
  }
  cat(str_flatten(unlist(list), collapse = "\n"))
}
#enumerkate(alphabet, n) # this one took me a long time

# introduction to random strings
string.probs <- function(string){
  gene <- strsplit(str_extract(string, "^\\w+(?=\\n)"), "")[[1]]
  probs <- as.numeric(str_split_1(str_extract(string, "(?<=\\n).*"), " "))
  
  vec <- rep(0, length(probs))
  for (i in 1:length(probs)) {
    vec[i] <- prod(case_match(gene, c("C", "G") ~    probs[i] /2, 
                              c("T", "A") ~ (1-probs[i])/2))
  }
  str_flatten(round(log10(vec), 4), collapse = " ")
}
#string.probs(string)

# HELPER FUNCTION: consolidate one pair of kmers
# takes seqs, breaks into starting and ending kmers, and does ONE ROUND of consolidating
# writes to global version of seqs, and if it consolidates, returns TRUE
consolidate.one.kmer <- function(sqs, k) {
  lis <- as.list(sqs)
  for (i in 1:length(sqs)) {
    lis[[i]] <- str_sub(sqs[i], c(1, str_length(sqs[i])-k+1),c(k, str_length(sqs[i])))
  }
  
  for (i in 1:length(lis)) {
    for (j in 1:length(lis)) {
      if (i == j) next
      if (lis[[i]][2] == lis[[j]][1]) {
        seqs[i] <<- sqs[i] <- str_remove(paste0(sqs[i], sqs[j]), lis[[i]][2])
        seqs <<- sqs[-j]
        return(TRUE)
      }
    }
  }
  return(FALSE)
}

# genome assembly as shortest superstring
assemble.reads <- function(fasta) {
  parse.fasta.v2(fasta)
  
  # determine a max and min k by finding half the length of the shortest read
  shortest <- seqs[which.min(str_length(seqs))]
  kmax <- str_length(shortest) - 1
  kmin <- round(str_length(shortest)/2, 0)
  klist <- kmin:kmax
  
  # loop through k's, doing as much consolidation at each k as possible
  for (k in klist) {
    if (length(seqs) == 1) break
    stillGoing <- TRUE
    while (stillGoing == TRUE) { # run consolidation fxn until it returns FALSE
      consolidate.one.kmer(seqs, k) -> stillGoing
    }
  }
  seqs
}
#assemble.reads(fasta)

# completing a tree
complete.tree <- function(n, string) {
  lis <- as.list(str_split(str_remove(string, "\\n$"), "\n")[[1]])
  tree <- list()
  
  # idea: put nodes into groups when connected by edges.
  # the number of new edges needed is the number of groups - 1, plus
  # the number of unconnected nodes
  for (i in 1:length(lis)) {
    vec <- as.numeric(str_split(lis[[i]], " ")[[1]])
    
    if (length(tree)==0) {tree[[1]] <- vec; next}
    
    pairInList <- FALSE   
    for (j in 1:length(tree)) { # scan existing groups for current nodes
      if (pairInList) break
      if (vec[1] %in% tree[[j]]) {
        for (k in 1:length(tree)) {
          if (j==k) next
          if (vec[2] %in% tree[[k]]) {
            tree[[j]] <- c(tree[[j]], tree[[k]])
            tree <- tree[-k]
            break
          }
        }
        tree[[j]] <- unique(c(tree[[j]], vec))
        pairInList <- TRUE
      } else if (vec[2] %in% tree[[j]]) { # repeat for second node in edge
        for (k in 1:length(tree)) {
          if (j==k) next
          if (vec[1] %in% tree[[k]]) {
            tree[[j]] <- c(tree[[j]], tree[[k]])
            tree <- tree[-k]
            break
          }
        }
        tree[[j]] <- unique(c(tree[[j]], vec)) 
        pairInList <- TRUE
      }
    }
    
    # if nodes not in list, create new slot for them
    if (!pairInList) {
      tree[[length(tree)+1]] <- vec
    }
  }
  loners <- sum(!(1:n %in% unlist(tree))) # find number of unconnected nodes
  length(tree) - 1 + loners
}
#complete.tree(n, string)

# partial permutations
partial.permutation <- function(n, k) return(prod(n:(n-k+1)) %% 1e6)
#partial.permutation(n, k)

# perfect matchings and RNA secondary structures
num.perfect.matchings <- function(fasta) {
  parse.fasta.v2(fasta)
  
  au <- str_count(seqs[1], "[AU]")/2
  cg <- str_count(seqs[1], "[CG]")/2
  
  format(factorialZ(au)*factorialZ(cg), scientific = F)
}
#num.perfect.matchings(fasta)

# enumerating oriented gene orderings
oriented.perms <- function(n) {
  total.possible <- prod(seq(2, (2*n), 2)) # n*(n-2)*(n-4)*...*2
  vec <- character(0)
  current <- numeric(n)
  cur.str <- character(1)
  
  while (length(vec) < total.possible) {
    current[] <- sample(1:n, n)
    current[] <- current*sample(c(-1, 1), n, replace = T)
    cur.str[] <- str_flatten(current, collapse = " ")
    if (!(cur.str %in% vec)) {
      vec <- c(vec, cur.str)
    } 
  }
  paste0(c(total.possible, vec), collapse = "\n")
  # note the output is too long to fit in the console, so write to a text file instead
}
#oriented.perms(n)

# longest increasing subsequence
subsequences <- function(permutation) {
  vec <- as.numeric(str_split(permutation, " ")[[1]])
  n <- length(vec)
  
  inc.subseqs <- vector("list", n) # despite the name, this creates a *list* of length n
  dec.subseqs <- vector("list", n)
  
  # increasing subsequences
  for (i in n:1) {
    best <- c(vec[i])
    # starting from the end, find the best increasing subsequence from that point
    # onward, then as you go backwards, append new values if they're smaller
    for (j in (i+1):n) {
      if (j>n) next
      if (vec[j] > vec[i] && length(inc.subseqs[[j]]) + 1 > length(best)) {
        best <- c(vec[i], inc.subseqs[[j]])
      }
    }
    inc.subseqs[[i]] <- best
  }
  best.inc <- str_flatten(inc.subseqs[[which.max(sapply(inc.subseqs, length))]], " ")
  
  # decreasing subsequences
  for (i in n:1) {
    best <- c(vec[i])
    
    for (j in (i+1):n) {
      if (j>n) next
      if (vec[j] < vec[i] && length(dec.subseqs[[j]]) + 1 > length(best)) {
        best <- c(vec[i], dec.subseqs[[j]])
      }
    }
    dec.subseqs[[i]] <- best
  }
  best.dec <- str_flatten(dec.subseqs[[which.max(sapply(dec.subseqs, length))]], " ")
  
  cat(c(best.inc, best.dec), sep = "\n")
}
#subsequences(permutation)
