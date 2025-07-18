library(tidyverse)

# complement a strand of DNA
revcomp <- function(string) {
  reversed <- strsplit(string, NULL)[[1]] |> rev() |> paste(collapse = "")
  strsplit(reversed, NULL)[[1]] |> 
    case_match("A" ~ "T",
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




