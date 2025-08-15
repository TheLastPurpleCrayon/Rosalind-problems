source("helper_functions.R")

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
#recc(n, k)

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
#compare_gc(string)

# counting point mutations
hamm <- function(string) {
  spl <- str_split(string, pattern = "\\n")[[1]]
  st1 <- spl[1]
  st2 <- spl[2]
  
  acc <- 0
  for (i in 1:str_length(st1)) acc <- acc + (str_sub(st1, i, i) != str_sub(st2, i, i))
  acc
}
#hamm(string)

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
#mendel(k, m, n)


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
#rna.to.prot(string)

# finding a motif in DNA
motif.pos <- function(string, motif) {
  locs <- str_locate_all(string, paste0("(?=(", motif, "))"))[[1]][,1]
  # have to use a lookahead to capture overlapping matches
  str_flatten(locs, collapse = " ")
}
#motif.pos(string, motif)