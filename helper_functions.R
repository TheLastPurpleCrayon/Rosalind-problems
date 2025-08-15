library(tidyverse)
library(gmp) # for big numbers
library(beepr)

# copied from Level3: complement a strand of DNA
revcomp <- function(string) {
  strsplit(string, NULL)[[1]] |>
    rev() |> 
    case_match("A" ~ "T", # could also do this with regex substitutions instead
               "T" ~ "A",
               "G" ~ "C",
               "C" ~ "G") |> 
    paste(collapse = "")
}

# HELPER FUNCTION: fasta parser
parse.fasta <- function(fasta) {
  spl <- str_split(fasta, pattern = ">")[[1]][-1]
  
  ids <<- str_extract(spl, pattern = "^.*(?=\\n)")
  seqs <<- str_extract(str_remove_all(spl, "\\n"), "[TACG]+")
  # after running, global environment should contain "ids" and "seqs" vars
}

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

# HELPER FUNCTION: tell me how long something takes to run and notify me when done
tictoc <- function(code) {
  tic <- Sys.time()
  cat(code)
  toc <- Sys.time()
  cat("\n")
  print(toc - tic)
  beep(3)
}

# HELPER FUNCTION: Hamming distance
hamming <- function(s1, s2) {
  v1 <- str_split_1(s1, "")
  v2 <- str_split_1(s2, "")
  acc <- 0
  for (i in 1:length(v1)) {
    acc <-  acc + (v1[i] != v2[i])
  }
  acc
}