revComp2 <- function(sequence) {
  model <- c("A", "a", "T", "t", "G", "g", "C", "c")
  names(model) <- c("T", "t", "A", "a", "C", "c", "G", "g")
  comp <- names(model)[match(strsplit(sequence, character())[[1]], model)]
  revComp <- rev(comp)
  paste(revComp, collapse = "")
}
