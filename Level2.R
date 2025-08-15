source("helper_functions.R")

# transcribing DNA into RNA
transcribe <- function(string) return(cat(str_replace_all(string, "T", "U")))
#transcribe(string)