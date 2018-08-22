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
    t(t(counts) / colSums(counts) * 10^6)
}

#' Rescale features
#'
#' Rescales each matrix row to the interval [0, 1].
#'
#' @name rescaleFeatures
#' @rdname rescaleFeatures
#' @aliases rescaleFeatures
#' @param counts matrix; a numeric matrix of counts.
#' @return The rescaled matrix.
#' @author Jason T. Serviss
NULL

#' @export

rescaleFeatures <- function(counts) {
  t(apply(counts, 1, function(x){
    (x - min(x)) / (max(x) - min(x))
  }))
}
