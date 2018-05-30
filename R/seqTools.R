#' \code{seqTools} package
#'
#' RNA sequencing tools.
#'
#'
#' @docType package
#' @name seqTools
#' @importFrom rlang .data
#' @importFrom dplyr "%>%"
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
