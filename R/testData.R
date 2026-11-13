#' Test counts data.
#'
#' @title Counts data used to demonstrate CIMseq package.
#' @docType data
#' @name testingCounts
#' @format matrix
#' \describe{
#'     \item{rownames}{Gene names}
#'     \item{colnames}{Samples. Singlets prefixed with "s" and multiplets "m".}
#' }
#' @usage testingCounts
#' @return Matrix of counts.
#' @examples
#' data(testingCounts)
NULL
"testingCounts"

#' Test counts meta data.
#'
#' @title Meta data for the testCounts dataset.
#' @docType data
#' @name testingMeta
#' @format tibble
#' \describe{
#'     \item{sample}{Sample ID}
#'     \item{cellNumber}{Indicates if samples are singlets or multiplets.}
#'     \item{cellTypes}{Indicates the cell types included in the sample.}
#' }
#' @usage testingMeta
#' @return Tibble with meta data.
#' @examples
#' data(testingMeta)
#'
NULL
"testingMeta"
