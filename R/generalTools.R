#' namedListToTibble
#'
#' Converts a named list to a long data frame.
#'
#' @name namedListToTibble
#' @rdname namedListToTibble
#' @author Jason T. Serviss
#' @param l List. The list to be converted.
#' @keywords namedListToTibble
#' @examples
#'
#' l <- list(a=LETTERS[1:10], b=letters[1:5])
#' output <- namedListToTibble(l)
#'
#' @export
#' @importFrom tibble tibble

namedListToTibble <- function(l) {
  #INPUT CHECKS
  if(!is.list(l)) {stop("This is not a list.")} #this won't catch a data frame
  if (length(names(l)) != length(l)) {
    stop("The list you submitted might not be named.")
  }

  if (!is.null(names(l[[1]]))) {
    #this handles when the list elements and the vector elements are both named
    #vector names and handling NA in vector names
    nas <- lapply(l, function(x) is.na(names(x)))
    ni <- gsub(".*\\.(.*)$", "\\1", names(unlist(l)))
    ni[unlist(nas)] <- NA
    #list element names
    n <- rep(names(l), lengths(l))
    tibble(names = n, inner.names = ni, variables = unname(unlist(l)))
  } else {
    #this handles when the list elements are named but the vector elements are unnamed
    n <- rep(names(l), lengths(l))
    tibble(names = n, variables = unname(unlist(l)))
  }
}

#' matrix_to_tibble
#'
#' Converts a matrix to a tibble without removing row names.
#'
#' @name matrix_to_tibble
#' @rdname matrix_to_tibble
#' @author Jason T. Serviss
#' @param data matrix; The matrix to be converted.
#' @param rowname character; Length 1 vector indicating the colname that
#'  rownames should have upon tibble conversion.
#' @param drop logical; indicates if rownames should be dropped.
#'  Default = FALSE.
#' @keywords matrix_to_tibble
#' @examples
#'
#' m <- matrix(rnorm(20), ncol = 2, dimnames = list(letters[1:10], LETTERS[1:2]))
#' output <- matrix_to_tibble(m)
#'
#' @export
#' @importFrom tibble as_tibble rownames_to_column

matrix_to_tibble <- function(data, rowname = "rowname", drop = FALSE) {
  if(!is.matrix(data)) stop("The 'data' argument is not a matrix")
  if(drop | is.null(rn)) as_tibble(data)

  data |>
  as.data.frame(stringsAsFactors = FALSE) |>
  add_column("{rowname}" := rownames(data), .before = 1) |>
  as_tibble(.name_repair = "universal_quiet")
}

#' tidySessionInfo
#'
#' there is a package for this now. pkg = sessioninfo, cmd = session_info()
#'
#' @name tidySessionInfo
#' @rdname tidySessionInfo
#' @author Jason T. Serviss
#' @examples
#' tidySessionInfo()
#'
#' @export
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom dplyr bind_rows bind_cols
#' @importFrom tidyr pivot_longer

#there is a package for this now. pkg = sessioninfo, cmd = session_info()
# tidySessionInfo <- function() {
#   l <- sessionInfo()
#   svs <- bind_cols(
#     R.version = l[[1]]$version.string,
#     platform = l[[1]]$platform,
#     locale = l[[3]],
#     running = l[[4]],
#     matprod = l[[8]],
#     BLAS = l[[9]],
#     LAPACK = l[[10]]
#   ) %>% gather(var, value) #change to pivot longer
#
#   basePkgs <- tibble(
#     var = rep("package_basePkgs", length(l[[5]])),
#     value = l[[5]]
#   )
#   otherPkgs <- tibble(
#     var = rep("package_otherPkgs", length(l[[6]])),
#     value = names(l[[6]])
#   )
#   loadedOnly <- tibble(
#     var = rep("package_loadedOnly"),
#     value = names(l[[7]])
#   )
#
#   bind_rows(svs, basePkgs, otherPkgs, loadedOnly)
# }
