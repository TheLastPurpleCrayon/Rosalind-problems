library(tidyverse)

# complement a strand of DNA
revcomp <- function(string) {
  strsplit(string, NULL)[[1]] |>
    rev() |> 
    case_match("A" ~ "T", # could also do this with regex substitutions instead
               "T" ~ "A",
               "G" ~ "C",
               "C" ~ "G") |> 
    paste(collapse = "")
}
revcomp(string)

# rabbits and recurrence relations
recc <- function(n, k) {
  numYoung <- 1
  numOld <- 0
  for (i in 1:(n-1)) {
    numOld <- numOld + numYoung
    numYoung <- (numOld-numYoung)*k
  }
  numOld + numYoung
}
recc(n, k)

# computing GC content
compare_gc <- function(string) {
  spl <- str_split(string, pattern = ">")[[1]]
  gc <- 0
  gc.best <- 0
  id.best <- ""
  gene.best <- ""
  for (i in 2:length(spl)) {
    id <- str_extract(spl[i], ".*(?=\\n)")
    gene <- str_replace_all(str_extract(spl[i], "(?<=\\n)(.|\\n)*"), "\\n", "")
    gc <- str_count(gene, pattern = "G|C")*100 / str_length(gene)
    if (gc > gc.best) {
      id.best <- id
      gene.best <- gene
      gc.best <- gc
    }
  }
  paste(id.best, gc.best, sep = "\n")
}
compare_gc(string)

# counting point mutations
hamm <- function(string) {
  spl <- str_split(string, pattern = "\\n")[[1]]
  st1 <- spl[1]
  st2 <- spl[2]
  
  acc <- 0
  for (i in 1:str_length(st1)) acc <- acc + (str_sub(st1, i, i) != str_sub(st2, i, i))
  acc
}
hamm(string)

# enumerating gene orders
perms <- function(n) {
  numPerms <- factorial(n)
  vec <- vector()
  
  while (length(vec) < numPerms) {
    cur <- sample(1:n, n, replace = F) |> as.character() |> str_flatten(collapse = " ")
    if (!(cur %in% vec)) vec <- append(vec, cur)
  }
  cat(paste0(c(numPerms, vec), collapse = "\n"))
}
perms(n)

# Mendel's first law
mendel <- function(k, m, n) {
  tot <- k + m + n
  homodomHomodom <- (k/tot)*(k-1)/(tot-1)
  homodomHet <- (k/tot)*(m/(tot-1))*2
  hetHet.75 <- (m/tot)*((m-1)/(tot-1))*0.75
  homodomHomorec <- (k/tot)*(n/(tot-1))*2
  hetHomorec.5 <- (m/tot)*(n/(tot-1))*2*0.5
  homodomHomodom + homodomHet + hetHet.75 + homodomHomorec + hetHomorec.5
}
mendel(k, m, n)


# translating RNA into protein
rna.to.prot <- function(string) {
  vec <- str_extract_all(string, "\\w{3}")[[1]]
  prot <- case_match(vec,
     # adapted from 'useful charts'/RNA_codon_chart.txt
     "UUU" ~ "F", "CUU" ~ "L", "AUU" ~ "I", "GUU" ~ "V",
     "UUC" ~ "F", "CUC" ~ "L", "AUC" ~ "I", "GUC" ~ "V",
     "UUA" ~ "L", "CUA" ~ "L", "AUA" ~ "I", "GUA" ~ "V",
     "UUG" ~ "L", "CUG" ~ "L", "AUG" ~ "M", "GUG" ~ "V",
     "UCU" ~ "S", "CCU" ~ "P", "ACU" ~ "T", "GCU" ~ "A",
     "UCC" ~ "S", "CCC" ~ "P", "ACC" ~ "T", "GCC" ~ "A",
     "UCA" ~ "S", "CCA" ~ "P", "ACA" ~ "T", "GCA" ~ "A",
     "UCG" ~ "S", "CCG" ~ "P", "ACG" ~ "T", "GCG" ~ "A",
     "UAU" ~ "Y", "CAU" ~ "H", "AAU" ~ "N", "GAU" ~ "D",
     "UAC" ~ "Y", "CAC" ~ "H", "AAC" ~ "N", "GAC" ~ "D",
     "UAA" ~ "",  "CAA" ~ "Q", "AAA" ~ "K", "GAA" ~ "E",
     "UAG" ~ "",  "CAG" ~ "Q", "AAG" ~ "K", "GAG" ~ "E",
     "UGU" ~ "C", "CGU" ~ "R", "AGU" ~ "S", "GGU" ~ "G",
     "UGC" ~ "C", "CGC" ~ "R", "AGC" ~ "S", "GGC" ~ "G",
     "UGA" ~ "",  "CGA" ~ "R", "AGA" ~ "R", "GGA" ~ "G",
     "UGG" ~ "W", "CGG" ~ "R", "AGG" ~ "R", "GGG" ~ "G"
  )
  str_flatten(prot)
}
rna.to.prot(string)


# inferring mRNA from protein
prot.to.rna <- function(string) {
  vec <- strsplit(string, "")[[1]]
  acc <- 1
  for (i in 1:length(vec)) {
    acc = acc * str_count(pattern = vec[i], 
      # string derived from flipping and concatenating the codon chart from before
      string = "FLIVFLIVLLIVLLMVSPTASPTASPTASPTAYHNDYHND_QKE_QKECRSGCRSG_RRGWRRG")
    if (acc > 1e6) acc <- acc %% 1e6 # numbers get too big for R otherwise
  }
  (acc*3) %% 1e6 # multiply by 3 to account for different stop codons
}
prot.to.rna(string)

# finding a motif in DNA
motif.pos <- function(string, motif) {
  locs <- str_locate_all(string, paste0("(?=(", motif, "))"))[[1]][,1]
    # have to use a lookahead to capture overlapping matches for some reason
  str_flatten(locs, collapse = " ")
}
motif.pos(string, motif)

# calculating expected offspring
expected <- function(string) {
  vec <- strsplit(string, " ")[[1]]
  j <- as.numeric(vec[1])
  k <- as.numeric(vec[2])
  l <- as.numeric(vec[3])
  m <- as.numeric(vec[4])
  n <- as.numeric(vec[5]) 
  o <- as.numeric(vec[6])
  (j*2) + (k*2) + (l*2) + (m*1.5) + (n*1) + (o*0)
}
expected(string)

# mortal fibonacci rabbits
library(gmp) # for big numbers
mort.fib <- function(n, m) {
  vec <- c(1, 0, rep(0, m-2)) # month 1
  vec <- as.bigz(vec)
  for (i in 2:n) {
    newAdults <- vec[1] # babies grow up
    vec[1] <- sum(vec)-vec[1] # adults reproduce
    for (j in length(vec):3) vec[j] <- vec[j-1] # adults move up the vector
    vec[2] <- newAdults # new adults are accounted for
  }
  sum(vec)
}
mort.fib(n, m)

# RNA splicing
rna.splice <- function(string) {
  spl <- str_split(string, pattern = ">")[[1]][-1]
  
  # splice, transcribe, translate
  temp <- str_split(spl[1], pattern = "\\n")[[1]]
  dna <- str_flatten(temp[-1])
  introns <- spl[2:length(spl)]
  for (i in 1:length(introns)) {
    intron <- str_split(introns[i], pattern = "\\n")[[1]][2]
    dna <- str_remove_all(dna, intron)
  }
  rna <- dna |> str_replace_all("T", "U")
  rna.to.prot(rna)
}
rna.splice(string)

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
subseq.idx(fasta)

# HELPER FUNCTION: fasta parser
parse.fasta <- function(fasta) {
  spl <- str_split(fasta, pattern = ">")[[1]][-1]
  
  ids <<- str_extract(spl, pattern = "^.*(?=\\n)")
  seqs <<- str_extract(str_remove_all(spl, "\\n"), "[TACG]+")
  # after running, global environment should contain "ids" and "seqs" vars
}

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
transition.transversion.ratio(fasta)

# finding a shared motif
comm.substring <- function(fasta) {
  parse.fasta(fasta)

  longest.match <- "_"
  for (k in 2:str_length(seqs[1])) {
    if ((k - str_length(longest.match)) >= 2) {
      #print("Done!")
      break # once longest substring is found, don't bother checking longer ones
    }
    #print(paste0("checking for ", k, "mers..."))
    kmers <- list()
    for (i in 1:length(seqs)) {
      vec <- str_sub(seqs[i], 1:(str_length(seqs[i])-(k-1)), k:str_length(seqs[i]))
      kmers <- c(kmers, list(vec))
    }
    for (i in 1:length(kmers[[1]])) {
      if (sum(str_detect(seqs, kmers[[1]][i])) == length(seqs)) {
        longest.match <- kmers[[1]][i]
        #print(paste0("new longest match: ", longest.match))
        break
      }
    }
  }
  longest.match
}
comm.substring(fasta)
# it turns out a better way to do this problem was to start from the shortest 
# sequence, test it, then test both its (k-1)mers, then its (k-2)mers, etc.
# that would probably speed things up a LOT compared to my bottom-up solution

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
enumerkate(alphabet, n) # this one took me a long time

# calculating protein mass
prot.mass <- function(aas) {
  spl <- strsplit(aas, "")[[1]]
  mass <- case_match(spl,
    # see protein_masses_(monoisotopic).txt, which
    "A" ~ "71.03711",
    "C" ~ "103.00919",
    "D" ~ "115.02694",
    "E" ~ "129.04259",
    "F" ~ "147.06841",
    "G" ~ "57.02146",
    "H" ~ "137.05891",
    "I" ~ "113.08406",
    "K" ~ "128.09496",
    "L" ~ "113.08406",
    "M" ~ "131.04049",
    "N" ~ "114.04293",
    "P" ~ "97.05276",
    "Q" ~ "128.05858",
    "R" ~ "156.10111",
    "S" ~ "87.03203",
    "T" ~ "101.04768",
    "V" ~ "99.06841",
    "W" ~ "186.07931",
    "Y" ~ "163.06333"
  )
  format(sum(as.numeric(mass)), nsmall = 5)
}
prot.mass(aas)

# consensus and profile
consensus.profile <- function(fasta) {
  parse.fasta(fasta)
  acgt <- c("A", "C", "G", "T")
  
  n <- str_length(seqs[1])
  prof <- matrix(0, nrow = 4, ncol = n, dimnames = list(acgt))
  for (i in 1:n) {              # each position
    for (j in 1:length(seqs)) { # each sequence
      for (k in 1:4) {          # each letter
        prof[k, i] <- prof[k, i] + (strsplit(seqs[j], "")[[1]][i] == acgt[k])
      }
    }
  }
  consensus <- character(n)
  for (l in 1:n) {
    consensus[l] <- acgt[which.max(prof[,l])]
  }
  cat(paste0(str_flatten(consensus), "\n", 
         "A: ", str_flatten(prof[1,], collapse = " "),
         "\nC: ", str_flatten(prof[2,], collapse = " "), 
         "\nG: ", str_flatten(prof[3,], collapse = " "), 
         "\nT: ", str_flatten(prof[4,], collapse = " ")))
}
consensus.profile(fasta)

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
string.probs(string)

# overlap graphs
adj.list.3 <- function(fasta) {
  parse.fasta(fasta)
  
  vec <- character(0)
  for (i in 1:length(ids)) {
    for (j in 1:length(ids)) {
      if (j == i) next
      if (str_extract(seqs[i], "\\w{3}$") == str_extract(seqs[j], "^\\w{3}")) {
        vec <- c(vec, paste0(ids[i], " ", ids[j]))
      }
    }
  }
  cat(str_flatten(vec, collapse = "\n"))
}
adj.list.3(fasta)

# HELPER FUNCTION: improved fasta parser
# this one works on protein fastas and doesn't get confused by ids containing [ACGT]
parse.fasta.v2 <- function(fa) {
  spl <- str_split(fa, pattern = ">")[[1]][-1]
  
  ids <<- str_extract(spl, pattern = "^.*(?=\\n)")
  for (i in 1:length(ids)) {
    spl[i] <- str_remove(spl[i], fixed(ids[i]))
    spl[i] <- str_remove_all(spl[i], "\\n")
  }
  seqs <<- spl
}

# finding a protein motif
motif.locations <- function(string) {
  prots <- strsplit(string, "\n")[[1]]
  
  output <- character(0)
  for (i in 1:length(prots)) {
    fasta <- read_file(paste0("http://www.uniprot.org/uniprot/", 
                              str_extract(prots[i], "^.{6}"), 
                              ".fasta"))
    parse.fasta.v2(fasta)
    
    locs <- str_flatten(str_locate_all(seqs[1], "(?=N[^P][ST][^P])")[[1]][,1], " ")
    if (locs != "") {
      output <- c(output, prots[i], locs)
    }
    
  }
  cat(str_flatten(output, "\n"))
}
motif.locations(string)

# HELPER FUNCTION: DNA trim
# takes a vector of DNA codons and returns those after START and before STOP, inclusive
dna.trim <- function(codons) {
  if (!("ATG" %in% codons)) return(character(0))
  while (codons[1] != "ATG") {
    codons <- codons[-1]
  }
  stops <- which(codons == "TAG" | codons == "TAA" | codons == "TGA")
  if (length(stops) == 0) return(character(0))
  valid <- str_flatten(codons[1:stops[1]])
  return(c(valid, dna.trim(codons[-1])))
}

# HELPER FUNCTION: DNA translator
# takes a DNA string and returns the translated peptide(s),
# respecting start and stop codons
dna.translate <- function(string) {
  codons <- str_extract_all(string, "\\w{3}")[[1]]
  if (length(codons) == 0) return(character(0))
  strings <- dna.trim(codons)
  if (length(strings) == 0) return(character(0))
  output <- character(0)
  for (i in 1:length(strings)){
    codons <- str_extract_all(strings[i], "\\w{3}")[[1]]
    aas <- case_match(codons, # see `useful charts`/DNA_codon_chart.txt
                      "TTT" ~ "F",   "CTT" ~ "L",   "ATT" ~ "I",   "GTT" ~ "V",
                      "TTC" ~ "F",   "CTC" ~ "L",   "ATC" ~ "I",   "GTC" ~ "V",
                      "TTA" ~ "L",   "CTA" ~ "L",   "ATA" ~ "I",   "GTA" ~ "V",
                      "TTG" ~ "L",   "CTG" ~ "L",   "ATG" ~ "M",   "GTG" ~ "V",
                      "TCT" ~ "S",   "CCT" ~ "P",   "ACT" ~ "T",   "GCT" ~ "A",
                      "TCC" ~ "S",   "CCC" ~ "P",   "ACC" ~ "T",   "GCC" ~ "A",
                      "TCA" ~ "S",   "CCA" ~ "P",   "ACA" ~ "T",   "GCA" ~ "A",
                      "TCG" ~ "S",   "CCG" ~ "P",   "ACG" ~ "T",   "GCG" ~ "A",
                      "TAT" ~ "Y",   "CAT" ~ "H",   "AAT" ~ "N",   "GAT" ~ "D",
                      "TAC" ~ "Y",   "CAC" ~ "H",   "AAC" ~ "N",   "GAC" ~ "D",
                      "TAA" ~ "",    "CAA" ~ "Q",   "AAA" ~ "K",   "GAA" ~ "E",
                      "TAG" ~ "",    "CAG" ~ "Q",   "AAG" ~ "K",   "GAG" ~ "E",
                      "TGT" ~ "C",   "CGT" ~ "R",   "AGT" ~ "S",   "GGT" ~ "G",
                      "TGC" ~ "C",   "CGC" ~ "R",   "AGC" ~ "S",   "GGC" ~ "G",
                      "TGA" ~ "",    "CGA" ~ "R",   "AGA" ~ "R",   "GGA" ~ "G",
                      "TGG" ~ "W",   "CGG" ~ "R",   "AGG" ~ "R",   "GGG" ~ "G")
    output <- c(output, str_flatten(aas))
  }
  output
}

# open reading frames
get.orfs <- function(fasta) {
  parse.fasta.v2(fasta)
  
  candidates <- character(6)
  for (i in 1:length(seqs)) {
    candidates[1] <- str_extract(seqs[i], "([ATCG]{3})*")
    candidates[2] <- str_extract(str_remove(seqs[i], "^[ATCG]{1}"), "([ATCG]{3})*")
    candidates[3] <- str_extract(str_remove(seqs[i], "^[ATCG]{2}"), "([ATCG]{3})*")
    rc <- revcomp(seqs[i])
    candidates[4] <- str_extract(rc, "([ATCG]{3})*")
    candidates[5] <- str_extract(str_remove(rc, "^[ATCG]{1}"), "([ATCG]{3})*")
    candidates[6] <- str_extract(str_remove(rc, "^[ATCG]{2}"), "([ATCG]{3})*")
  }
  
  orfs <- character(0)
  for (i in 1:length(candidates)) {
    orfs <- c(orfs, dna.translate(candidates[i]))
  }
  cat(paste0(unique(orfs), collapse = "\n"))
}
get.orfs(fasta)

# HELPER FUNCTION: complement
# return the base which matches the given DNA base (so no U)
complement <- function(base) {
  if (base[1] == "A") return("T")
  else if (base[1] == "T") return("A")
  else if (base[1] == "C") return("G")
  else if (base[1] == "G") return("C")
  else return(NA)
}

# locating restriction sites 
rev.palindromes <- function(fasta) {
  parse.fasta.v2(fasta)
  vec <- strsplit(seqs[1], "")[[1]]
  # filter: at each position, check if its complement is 4-12 bases away
  strings <- character(0)
  indices <- character(0)
  for (i in 1:length(vec)) {
    for (j in 3:11) {
      if ((i+j) > length(vec)) next
      if (complement(vec[i]) == vec[i+j]) {
        strings <- c(strings, str_flatten(vec[i:(i+j)]))
        indices <- c(indices, paste0(c(i, j+1), collapse = " "))
      }
    }
  }
  
  # find reverse palindromes from pre-filtered sequences
  output <- character(0)
  for (i in 1:length(strings)) {
    if (strings[i] == revcomp(strings[i])) output <- c(output, indices[i])
  }
  cat(paste0(unique(output), collapse = "\n"))
}
rev.palindromes(fasta)

# independent alleles
doub.heterozygous <- function(k, n) return(1-pbinom(n-1, 2^k, 1/4))
# crossing *anything* with a double heterozygote yields a double heterozygote 
# exactly 1/4 of the time no matter what the genotype of the other parent
doub.heterozygous(k, n)


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
assemble.reads(fasta)

# # longest increasing subsequence
# find.subsequences <- function(permutation) {
#   vec <- as.numeric(str_split(permutation, " ")[[1]])
#   n <- length(vec)
#   
#   # find increasing subsequence
#   longest <- numeric(0)
#   for (i in 1:(n-1)) {
#     working <- vec[i]
#     for (j in (i+1):n) {
#       if (vec[j] > working[length(working)]) {working <- c(working, vec[j])}
#     }
#     if (length(working) > length(longest)) {longest <- working}
#   }
#   increasing <- str_flatten(longest, " ")
#   
#   # find decreasing subsequence
#   longest <- numeric(0)
#   for (i in 1:(n-1)) {
#     working <- vec[i]
#     for (j in (i+1):n) {
#       if (vec[j] < working[length(working)]) {working <- c(working, vec[j])}
#     }
#     if (length(working) > length(longest)) {longest <- working}
#   }
#   decreasing <- str_flatten(longest, " ")
#   
#   cat(paste0(c(increasing, decreasing), collapse = "\n"))
# }
# find.subsequences(permutation)
# 
# # helper function
# # take the first item in the given vector
# #   remove all items in the vector smaller than that item
# #   for all remaining items,
# #     take the first item, and remove all values smaller than it
# #     but also take the second item, and remove all values smaller than it
# 
# remove.and.crunch <- function(string, i) {
#   vec <- str_split(string, " ")[[1]]
#   if (length(vec) == 1) return(vec)
#   ith <- vec[i]
#   for (j in 1:i) {
#     vec <- vec[j]
#   }
#   vec <- vec[which(vec >= ith)]
#   return(c(vec[i], remove.and.crunch(str_flatten(vec, " "), i)))
# }
# remove.and.crunch(vec, 1)
# 
# # for each value
# # compile a list of increasing substrings
# #   for j in i:n
# #     
# # take the longest
# 
# fxn <- function(vec, i) {
#   if (length(vec) <= 1) return(vec)
#   ith <- vec[i]
#   remaining <- vec[-(1:i)]
#   reamining <- remaining[which(remaining >= ith)]
#   for (j in 1:length(remaining)) {
#     res <- c(ith, fxn(remaining, j))
#     return(res)
#   }
# }
# 
# fxn.wrapper <- function(vec, max.i) {
#   accumulator <- character(0)
#   for (k in 1:max.i) {
#     accumulator <- list(accumulator, fxn(vec, k))
#   }
#   accumulator
# }
# fxn.wrapper(vec, 4)
# 
# 
# 
# 
# # HELPER FUNCTIONS: randomize an increasing/decreasing subsequence
# find.increasing.subsequence <- function(vec, numIter) {
#   n <- length(vec)
#   longest <- numeric(0)
#   
#   for (i in 2:n) {
#     # take numIter tries to randomize an increasing subsequence of size i 2:n
#     samp <- numeric(i)
#     for (b in 1:numIter) {
#       samp[] <- sample(vec, size = i)
#       samp[] <- vec[which(vec %in% samp)] # this puts the items in the original order
#       if (all(diff(samp) > 0)) { # if sample is strictly increasing, store it
#         longest <- samp
#         break
#       } else if (b == numIter) { # if the for loop completes with no match, be done
#         return(longest)
#       }
#     }
#   }
#   stop("No increasing subsequence found")
# }
# find.decreasing.subsequence <- function(vec, numIter) {
#   n <- length(vec)
#   longest <- numeric(0)
#   
#   for (i in 2:n) {
#     # take numIter tries to randomize an increasing subsequence of size i 2:n
#     samp <- numeric(i)
#     for (b in 1:numIter) {
#       samp[] <- sample(vec, size = i)
#       samp[] <- vec[which(vec %in% samp)] # this puts the items in the original order
#       if (all(diff(samp) < 0)) { # if sample is strictly decreasing, store it
#         longest <- samp
#         break
#       } else if (b == numIter) { # if the for loop completes with no match, be done
#         return(longest)
#       }
#     }
#   }
#   stop("No decreasing subsequence found")
# }
# 
# generate.subsequences <- function(string, numIter = 1e6) {
#   string <- str_remove_all(string, "\\r") # manage weird whitespace from file read
#   string <- str_remove(string, "^.*\\n") # take off the first line (permutation length)
#   
#   vec <- as.numeric(str_split(string, " ")[[1]])
#   
#   inc <- find.increasing.subsequence(vec, numIter) |> str_flatten(collapse = " ")
#   dec <- find.decreasing.subsequence(vec, numIter) |> str_flatten(collapse = " ")
#   
#   cat(paste0(c(inc, dec), collapse = "\n"))
# }
# Sys.time(); generate.subsequences(string, numIter = 1e6); Sys.time()


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
complete.tree(n, string)


# partial permutations
partial.permutation <- function(n, k) return(prod(n:(n-k+1)) %% 1e6)
partial.permutation(n, k)

# perfect matchings and RNA secondary structures
num.perfect.matchings <- function(fasta) {
  parse.fasta.v2(fasta)
  
  au <- str_count(seqs[1], "[AU]")/2
  cg <- str_count(seqs[1], "[CG]")/2
  
  format(factorialZ(au)*factorialZ(cg), scientific = F)
}
num.perfect.matchings(fasta)

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
oriented.perms(n)

# counting subsets
num.subsets <- function(n) as.bigz(2)^n %% 1e6
num.subsets(n)
