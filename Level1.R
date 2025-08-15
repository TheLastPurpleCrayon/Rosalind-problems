source("helper_functions.R")

# counting DNA nucleotides
count.nucleotides <- function(string) {
  As <- str_count(string, "A")
  Cs <- str_count(string, "C")
  Gs <- str_count(string, "G")
  Ts <- str_count(string, "T")
  cat(As, Cs, Gs, Ts, sep = " ")
}
#count.nucleotides(string)
