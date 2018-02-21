
revComp <- function(sequence) {
  as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(sequence)))
}

#test
sequence <- c("A", "T", "G", "C")
revComp(sequence)

#probably faster with
revComp2 <- function(sequence) {
  dplyr::case_when(
    sequence == "A" ~ "T",
    sequence == "T" ~ "A",
    sequence == "G" ~ "C",
    sequence == "C" ~ "G"
  )
}

revComp2(sequence)
