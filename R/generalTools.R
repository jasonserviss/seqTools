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
#' @importFrom tibble as_tibble add_column
#' @importFrom rlang :=

matrix_to_tibble <- function(data, rowname = "rowname", drop = FALSE) {
  if(!is.matrix(data)) stop("The 'data' argument is not a matrix")
  if(drop) as_tibble(data)

  data |>
  as.data.frame(stringsAsFactors = FALSE) |>
  add_column("{rowname}" := rownames(data), .before = 1) |>
  as_tibble(.name_repair = "universal_quiet")
}

#' divide_by
#'
#' Facilitates division in pipes and "avoids wrong number of arguments to"
#' complaint by check. Adopted from magrittr:
#' \url{https://github.com/tidyverse/magrittr/blob/master/R/aliases.R#L93-L96}
#'
#' @name divide_by
#' @rdname divide_by
#' @author Jason T. Serviss
#' @export

divide_by <- `/`

#' multiply_by
#'
#' Facilitates multiplication in pipes and "avoids wrong number of arguments to"
#' complaint by check. Adopted from magrittr:
#' \url{https://github.com/tidyverse/magrittr/blob/master/R/aliases.R#L82-L85}
#'
#' @name multiply_by
#' @rdname multiply_by
#' @author Jason T. Serviss
#' @export

multiply_by <- `*`

#' col40
#'
#' Diverging color palette with 40 colors.
#'
#' @name col40
#' @rdname col40
#' @author Jason T. Serviss
#' @keywords col40
#' @examples
#'
#' col40()
#'
#' @export
NULL

col40 <- function() {
  c(
  "#1c54a8", "#f63b32", "#00C2A0", "#FFDBE5", "#e0c48c", "#63FFAC", "#663000",
  "#e0a81c", "#385438", "#609060", "#6A3D9A", "#548495", "#A30059", "#8FB0FF",
  "#997D87", "#4FC601", "#8ca8c4", "#3B5DFF", "#BA0900", "#DDEFFF", "#7B4F4B",
  "#A1C299", "#0AA6D8", "#00846F", "#FFB500", "#C2FFED", "#A079BF", "#C0B9B2",
  "#C2FF99", "#00489C", "#6F0062", "#EEC3FF", "#922329", "#FFF69F", "#FF8A9A",
  "#B05B6F", "#7900D7", "#BC23FF", "#9B9700", "#0089A3"
  )
}
