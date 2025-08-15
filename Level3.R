source("helper_functions.R")

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
#revcomp(string)