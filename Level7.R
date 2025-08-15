source("helper_functions.R")

# counting subsets
num.subsets <- function(n) as.bigz(2)^n %% 1e6
#num.subsets(n)

# counting phylogenetic ancestors
ancestor.count <- function(n) return(n-2)
#ancestor.count(n)

# speeding up motif finding
kmp.failure.array <- function(fasta) {
  parse.fasta.v2(fasta)
  string <- str_split(seqs[1], "")[[1]]
  
  vec <- numeric(length(string))
  for (k in 2:length(vec)) {
    # if the previous entry is 0, only check if current should be 1
    # if the previous entry is nonzero, check values from 1+previous value downward
    if (vec[k-1] == 0) {
      vec[k] <- string[k] == string[1]
    } else {
      for (i in (vec[k-1]+1):1) {
        if (str_flatten(string[(k-i+1):k]) == str_flatten(string[1:i])) {
          vec[k] <- i
          break
        }
      }
    }
  }
  str_flatten(vec, collapse = " ")
}
#kmp.failure.array(fasta)

# matching random motifs
random.motif.prob <- function(n, gc, string) {
  vec <- str_split_1(string, "")
  p <- 1
  
  for (i in 1:length(vec)) {
    if (vec[i] %in% c("G", "C")) p <- p * (gc/2)
    else p <- p * ((1-gc)/2)
  }
  
  1-dbinom(0, n, prob = p)
}
#random.motif.prob(n, gc, string)

# HELPER FUNCTION: given an alphabet, convert a string into its ABCD equiv for sorting
my.match <- function(alphabet, string) {
  s <- str_split_1(string, "")
  vec <- character(length(s))
  for (i in 1:length(s)) {
    vec[i] <- LETTERS[which(alphabet == s[i])]
  }
  str_flatten(vec, "")
}

# ordering of strings of varying length lexicographically
enumerate.in.order <- function(string, n) {
  alphabet <- str_split_1(string, " ")
  entries <- character(0)
  
  for (i in 1:n) {
    vec <- character(0)
    max <- length(alphabet)^i
    while (length(vec) < max) { # there's probably a faster way than randomization
      candidate <- str_flatten(sample(alphabet, size = i, replace = T), "")
      if (!(candidate %in% vec)) vec <- c(vec, candidate)
    }
    entries <- c(entries, vec)
  }
  
  # sorting requires a conversion to the normal alphabet
  numbers <- character(length(entries))
  for (i in 1:length(entries)) {
    numbers[i] <- my.match(alphabet, entries[i])
  }
  tib <- tibble(entries = entries, numbers = numbers) # work around lack of dicts
  output <- tib |> 
    arrange(numbers) |> 
    pull(entries)
  
  cat(output, sep = "\n", file = file)
}
#enumerate.in.order(string, n)

# HELPER FUNCTION: make a reversal of a permutation
arb.rev <- function(perm, position, halfwidth) {
  # perm: the permutation to make a reversal in
  # position: the middle of the reversal. supports both integers (n) and n.5
  # halfwidth: the full width of the reversal divided by 2
  bounds <- position + c(-1, 1)*halfwidth
  start <- ceiling(bounds[1])
  end <- floor(bounds[2])
  perm[start:end] <- perm[end:start]
  perm
}

# HELPER FUNCTION: given two vectors of the same numbers, return the total
# distance between the two, measured as the sum of the numbers of spaces
# each element had to move to get from one to the other
totalDistance <- function(v1, v2) {
  n <- length(v1)
  dist <- 0
  for (i in 1:n) {
    dist <- dist + abs(which(v2 == v1[i]) - i)
  }
  dist
}

# reversal distance
reversal.distance <- function(string) {
  # more parsing than usual for this problem
  vec <- str_split_1(string, "\n")
  origs <- character(0)
  revds <- character(0)
  for (i in 1:length(vec)) {
    if (i %% 3 == 1) {
      origs <- c(origs, vec[i])
    } else if (i %% 3 == 2) {
      revds <- c(revds, vec[i])
    }
  }
  
  output <- numeric(0)
  for (i in 1:length(origs)) {
    orig <- as.numeric(str_split_1(origs[i], " "))
    revd <- as.numeric(str_split_1(revds[i], " "))
    n <- length(orig)
    
    totDist <- totalDistance(orig, revd)
    
    
    old <- list(revd)
    cand <- numeric(n)
    
    distances <- 50
    for (b in 1:3e3){
      numIters <- 0
      numRevs <- 0
      numIterationsSinceImprovement <- 0
      totDist <- totalDistance(orig, revd)
      working <- revd
      while(totDist > 0 && numRevs <= min(distances)) {
        
        # randomly do a reversal, taking pains to make sure it's valid
        pos <- sample(seq(1.5, n-0.5, by=0.5), 1)
        maxhw <- min(floor(pos-0.5), floor(n+0.5-pos))
        hw <- sample(1:maxhw, 1)
        cand <- arb.rev(working, position = pos, halfwidth = hw)
        
        # if the reversal moves us closer to the original, go from there, discard o/w
        newDist <- totalDistance(orig, cand)
        if (newDist < totDist) {
          old[[length(old)+1]] <- working
          working <- cand
          totDist <- newDist
          numRevs <- numRevs + 1
          numIterationsSinceImprovement <- 0
        } else {
          numIterationsSinceImprovement <- numIterationsSinceImprovement + 1
        }
        
        # failsafe: if no improvement in last 100 iterations, go back one reversal
        if (numIterationsSinceImprovement > 500) {
          if (length(old)>0) working <- old[[length(old)]]
          old <- old[-length(old)]
          numIterationsSinceImprovement <- 0
        }
        # failsafe 2: if no movement in 10000 iterations, skip and move on
        numIters <- numIters + 1
        if (numIters > 1e4) {
          numRevs <- Inf
          break
        }
        #cat(i, totDist, numRevs, min(distances), "\n")
      }
      distances <- c(distances, numRevs)
    }
    # as.numeric(names(table(distances)[which.max(table(distances))]))
    output <- c(output, min(distances))
  }
  cat(output, " ")
  # sometimes this function does not find the minimum number of reversals, sometimes it
  # finds a number *lower* than the minimum somehow. but at least it doesn't
  # require 20 mins of enumeration like the intended solution does!
} 
#reversal.distance(string)

# maximum matchings and RNA secondary structures
maximum.matchings <- function(fasta) {
  parse.fasta.v2(fasta)
  string <- seqs[1]
  
  cg <- c(str_count(string, "C"), str_count(string, "G"))
  au <- c(str_count(string, "A"), str_count(string, "U"))
  
  max.cg <- max(cg)
  min.cg <- min(cg)
  out <- factorialZ(max.cg)/factorialZ(max.cg - min.cg)
  max.au <- max(au)
  min.au <- min(au)
  out <- out*(factorialZ(max(au))/factorialZ(max.au - min.au))
  
  out
}
#maximum.matchings(fasta)

# HELPER FUNCTION: calculate the Catalan numbers
catalans <- function(num) {
  nums <- c("0" = 1, "1" = 1)
  for (n in 2:num) {
    acc <- 0
    for (k in 1:n) {
      acc <- acc + (nums[as.character(k-1)] * nums[as.character(n-k)])
    }
    nums[as.character(n)] <- acc
  }
  nums
}

# HELPER FUNCTION: recursion for counting the number of noncrossing matches
count.matches <- function(s) {
  #print(vec)
  if (nchar(s) == 0) return(1)
  if (!is.null(memo[[s]])) return(memo[[s]])
  
  poss.matches <- str_locate_all(s, matches[substr(s, 1, 1)])[[1]][,1]
  poss.matches <- poss.matches[which(poss.matches %% 2 == 0)]
  
  if (length(poss.matches) == 0) return(0)
  
  acc <- 0
  for (i in 1:length(poss.matches)) {
    if (poss.matches[i] == 2) {
      betweens <- ""
    } else {
      betweens <- substr(s, 2, poss.matches[i] - 1)
    }
    
    if (poss.matches[i] == nchar(s)) {
      outsides <- ""
    } else {
      outsides <- substr(s, poss.matches[i] + 1, nchar(s))
    }
    
    count.betw <- count.matches(betweens)
    count.outs <- count.matches(outsides)
    
    acc <- acc + (count.betw*count.outs) %% 1e6
  }
  memo[[s]] <<- acc
  acc
}

# catalan numbers and RNA secondary structures
noncrossing.matchings <- function(fasta) {
  parse.fasta.v2(fasta)
  string <- seqs[1]
  
  matches <<- c("A" = "U", "U" = "A", "C" = "G", "G" = "C")
  memo <<- list()
  count.matches(string) %% 1e6
}
#noncrossing.matchings(fasta)

# creating a distance matrix
dist.mat <- function(fasta) {
  parse.fasta.v2(fasta)
  n <- length(seqs)
  mat <- matrix(NA, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i==j) {
        mat[i, j] = 0
      } else {
        mat[i, j] = sum(str_split_1(seqs[i], "") != str_split_1(seqs[j], ""))
      }
    }
  }
  mat <- mat/nchar(seqs[1])
  prmatrix(mat, rowlab = rep("", n), collab = rep("", n))
}
#dist.mat(fasta)

# HELPER FUNCTION: modify enumerkate() above to return its output rather than
# printing it
enumerkate.return <- function(alphabet, n) {
  alph <- str_remove_all(alphabet, "\\s")
  
  list <- list(strsplit(alph, "")[[1]])
  for (i in 1:(n-1)){
    for (j in length(list[[1]]):1) {
      list[[j]] <- paste0(list[[1]][j], strsplit(alph, "")[[1]])
    }
    list <- list(unlist(list))
  }
  unlist(list)
}

# k-mer composition
tetranucleotide.composition <- function(fasta) {
  parse.fasta.v2(fasta)
  fourmers <- enumerkate.return("ACGT", 4)
  
  pat <- paste0("(?=(", fourmers, "))")
  
  #out <- str_count(seqs[1], fourmers) 
  # for some reason, str_count doesn't count overlapping matches correctly
  out <- sapply(str_locate_all(seqs[1], pat), length)/2
  
  cat(out, collapse = " ")
}
#tetranucleotide.composition(fasta)

# error correction in reads
corrections <- function(fasta) {
  parse.fasta.v2(fasta)
  
  seen.once <- character(0)
  seen.twice <- character(0)
  for (i in 1:length(seqs)) {
    if (seqs[i] %in% seen.once || revcomp(seqs[i]) %in% seen.once) {
      seen.twice[length(seen.twice)+1] <- seqs[i]
      toRemove <- which(seen.once == seqs[i])
      if (length(toRemove) == 0) toRemove <- which(seen.once == revcomp(seqs[i]))
      seen.once <- seen.once[-toRemove]
    } else if (seqs[i] %in% seen.twice || revcomp(seqs[i]) %in% seen.twice) {
      # do nothing
    } else {
      seen.once[length(seen.once)+1] <- seqs[i]
    }
  }
  
  o <- character(0)
  hamming.vec <- Vectorize(hamming)
  for (j in 1:length(seen.once)) {
    ref <- seen.twice[which(hamming.vec(seen.once[j], seen.twice) == 1)]
    if (length(ref) == 0) {
      ref <- seen.twice[which(hamming.vec(revcomp(seen.once[j]), seen.twice) == 1)]
    }
    o <- c(o, paste0(seen.once[j], "->", ref))
  }
  cat(o, sep = "\n")
}
#corrections(fasta)

# HELPER FUNCTION: recursive search for subsequences
# Uses memoization to speed up a recursion that ultimately builds back to front
comm.subseq.search <- function(idxOfLast1, idxOfLast2) {
  
  # check if problem has been done already
  key <- paste(idxOfLast1, idxOfLast2, sep = ",")
  if (exists(key, envir = memo)) {
    return(get(key, envir = memo))
  } 
  
  best_subseq = ""
  
  for (i in 1:length(acgt)) {
    s1match <- regexpr(acgt[i], substr(seqs[1], idxOfLast1+1, nchar(seqs[1])))[[1]]
    s2match <- regexpr(acgt[i], substr(seqs[2], idxOfLast2+1, nchar(seqs[2])))[[1]]
    if ((s1match != -1) && (s2match != -1)) {
      cand <- paste0(acgt[i], 
                     comm.subseq.search(idxOfLast1+s1match, idxOfLast2+s2match))
      if (nchar(cand) > nchar(best_subseq)) best_subseq <- cand
    }
  }
  
  assign(key, best_subseq, envir = memo)
  return(best_subseq)
}

# finding a shared spliced motif
longest.common.subsequence <- function(fasta) {
  parse.fasta.v2(fasta)
  
  # going to use base R string manipulation as a challenge
  # str_length() -> nchar()
  # str_split_1 and vec[i] -> substr(string, i, i)
  
  acgt <<- c("A", "C", "G", "T")
  
  cat(comm.subseq.search(0, 0))
}
#memo <- new.env(hash = T, parent = emptyenv())
#longest.common.subsequence(fasta)
#tictoc(longest.common.subsequence(fasta))

