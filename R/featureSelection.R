#' nTopVar
#'
#' Facilitates gene selection using variance.
#'
#' Returns the index for the n genes (rows) with the maximum
#' variance in the counts object. Counts per million (CPM) should be used for
#' the calculation.
#'
#' @name nTopVar
#' @rdname nTopVar
#' @aliases nTopVar
#' @param cpm matrix; Matrix containing cpm values.
#' @param n Number of genes to select.
#' @return A numeric vector containing the indices of selected genes.
#' @author Jason T. Serviss
#' @keywords nTopVar
#'
NULL

#' @rdname nTopVar
#' @importFrom stats var
#' @importFrom matrixStats rowVars
#' @export

nTopVar <- function(cpm, n) {
  n <- min(n, dim(cpm)[1])
  rv = rowVars(cpm)
  order(rv, decreasing = TRUE)[1:n]
}

#' nTopMax
#'
#' Facilitates gene selection prior to unsupervised clustering.
#'
#' Returns the index for the n genes (rows) with the maximum
#' expression in the spCounts object. The expression matrix in
#' the counts.cpm slot is used for the calculation.
#'
#' @name nTopMax
#' @rdname nTopMax
#' @aliases nTopMax
#' @param cpm matrix; Matrix containing cpm values.
#' @param n Number of genes to select.
#' @return A numeric vector containing the indices of selected genes.
#' @author Jason T. Serviss
#' @keywords nTopMax
#'
NULL

#' @rdname nTopMax
#' @export

nTopMax <- function(cpm, n) {
  n <- min(n, dim(cpm)[1])
  rv <- apply(cpm, 1, max)
  select <- order(rv, decreasing = TRUE)[1:n]
  return(select)
}
