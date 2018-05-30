
revComp <- function(sequence) {
  as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(sequence)))
}

#test
sequence <- c("A", "T", "G", "C")
revComp(sequence)

#even faster with
revComp2 <- function(sequence) {
  model <- c("A", "a", "T", "t", "G", "g", "C", "c")
  names(model) <- c("T", "t", "A", "a", "C", "c", "G", "g")
  comp <- names(model)[match(strsplit(sequence, character())[[1]], model)]
  revComp <- rev(comp)
  paste(revComp, collapse = "")
}

revComp2(sequence)
