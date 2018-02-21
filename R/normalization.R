#' cpm
#'
#' Calculates counts per million (cpm) using a gene expression counts matrix
#' as input.
#'
#' @name cpm
#' @rdname cpm
#' @aliases cpm
#' @param counts matrix; a numeric matrix of counts.
#' @return A matrix of cpm values.
#' @author Jason T. Serviss
NULL

#' @export

cpm <- function(counts) {
    norm.fact <- colSums(counts)
    counts.cpm <- t(apply(counts, 1, .norm, n = norm.fact))
}

.norm <- function(x, n) {
    x / n * 1000000 + 1
}
