#' cpm
#'
#' Calculates counts per million (cpm) using a scRNAseq counts matrix
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
    #do I want input checks here?
    t(t(counts) / colSums(counts) * 10^6)
}

#' cpm.log2
#'
#' Calculates counts per million (cpm) using a scRNAseq counts matrix
#' as input.
#'
#' @name cpm.log2
#' @rdname cpm.log2
#' @aliases cpm.log2
#' @param counts matrix; a numeric matrix of counts.
#' @return A matrix of cpm values.
#' @author Jason T. Serviss
NULL

#' @export

cpm.log2 <- function(counts) {
  log2(cpm(counts) + 1)
}


#' Rescale features
#'
#' Rescales each matrix row to the interval \\[0, 1\\].
#'
#' @name rescale.features
#' @rdname rescale.features
#' @aliases rescale.features
#' @param counts matrix; a numeric matrix of scRNAseq counts.
#' @return The rescaled matrix.
#' @author Jason T. Serviss
NULL

#' @export

rescale.features <- function(counts) {
  t(apply(counts, 1, function(x){
    (x - min(x)) / (max(x) - min(x))
  }))
}

#' Pearsons distance
#'
#'
#' Calculates 1-Pearson's correlation as a distance.
#'
#' @name pearsons.dist
#' @rdname pearsons.dist
#' @aliases pearsons.dist
#' @param cpm matrix; Counts per million.
#' @param select A numeric vector indicating the indexes of genes to include.
#' @return A matrix containing the mean value for each gene for each
#'    classification group.
#' @author Jason T. Serviss
#' @keywords pearsons.dist
NULL

#' @rdname pearsons.dist
#' @importFrom stats cor as.dist
#' @export

pearsons.dist <- function(cpm, select) {
  as.dist(1 - cor(cpm[select, ], method = "p"))
}
