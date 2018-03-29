
#' my.ks.test
#'
#' A faster version of ks.test when ties are present. Runs ks.test with args
#' exact = TRUE and alternative = "two.sided". Errors if NAs or non-numeric data
#' are present in x or y. In addition, errors if length(x) or length(y) < 1.
#' Note that this function calls C code in the base R package which isn't
#' typically "allowed" and, in addition, it is specifically stated that
#' "Packages should not make .C/.Call/.External/.Fortran calls to a base
#' package. They are not part of the API, for use only by R itself and subject
#' to change without notice."
#'
#' @name my.ks.test
#' @rdname my.ks.test
#' @aliases my.ks.test
#' @param x Numeric; a numeric vector of data values.
#' @param y Numeric; a numeric vector of data values.
#' @param check Logical; should input checks be run.
#' @param ... Arguments to pass on.
#' @return A list with the statistic as element 1 and the p.value as element 2.
#' @author Jason T. Serviss
#' @examples
#' my.ks.test(rnorm(100), rnorm(100))
#' my.ks.test(1:100, 101:200)
#' my.ks.test(rnorm(100), rnorm(50))
#' my.ks.test(1:100, 51:150)
#' \dontrun{
#'  x <- rnorm(1000)
#'  y <- rnorm(1000)
#'  microbenchmark(microbenchmark(ks.test(x, y), my.ks.test(x, y)))
#'  x <- 1:1000
#'  y <- 501:1500
#'  microbenchmark(microbenchmark(ks.test(x, y), my.ks.test(x, y)))
#' }
#' @export

my.ks.test <- function (x, y, check = FALSE, ...){
  
  if(check) {
    if(any(is.na(x))) {stop("NAs in x")}
    if(any(is.na(y))) {stop("NAs in y")}
    if(any(!is.numeric(x))) {stop("Non-numeric values in x")}
    if(any(!is.numeric(y))) {stop("Non-numeric values in y")}
    if(length(x) < 1L) {stop("Not enough x data")}
    if(length(y) < 1L) {stop("Not enough y data")}
  }
  
  PVAL <- NULL
  TIES <- FALSE
  
  n <- length(x)
  n.x <- as.double(n)
  n.y <- length(y)
  n <- n.x * n.y / (n.x + n.y)
  w <- c(x, y)
  z <- cumsum(ifelse(order(w) <= n.x, 1 / n.x, -1 / n.y))
  
  if (length(unique(w)) < (n.x + n.y)) {
    TIES <- TRUE
    z <- z[c(which(diff(sort(w)) != 0), n.x + n.y)]
  }
  
  STATISTIC <- max(abs(z))
  
  if (!TIES) {
    PVAL <- 1 - .Call(stats:::C_pSmirnov2x, STATISTIC, n.x, n.y)
  }
  
  if (is.null(PVAL)) {
    pkstwo <- function(x, tol = 1e-06) {
      x <- as.double(x)
      p <- rep(0, length(x))
      IND <- which(x > 0)
      if (length(IND))
      p[IND] <- .Call(stats:::C_pKS2, p = x[IND], tol)
      p
    }
    PVAL <- 1 - pkstwo(sqrt(n) * STATISTIC)
  }
  PVAL <- min(1, max(0, PVAL))
  return(list(STATISTIC, PVAL))
}

#' KStest
#'
#' Runs my.ks.test (giving equivilent results to ks.test) comparing every gene
#' in exp, for all possible combinations of classes.
#'
#' @name KStest
#' @rdname KStest
#' @aliases KStest
#' @param counts Matrix; counts matrix with samples as columns and genes as rows.
#'  Should have both colnames and rownames.
#' @param classes Character; vector with length ncol(exp) indicating the class
#'  of each sample.
#' @param cores Numeric; indicates the number of cores to run the analysis on.
#' @return A tibble with one row per gene and class combination and the
#'  corresponding test statistic and p.value.
#' @author Jason T. Serviss
#' @examples
#' #nothing here yet
#' @export
#' @importFrom parallel mclapply
#' @importFrom tibble as_tibble
#' @importFrom dplyr bind_rows
#' @importFrom magrittr "%>%"
#' @importFrom purrr map2_dfr
#' @importFrom stats setNames
#' @importFrom utils combn

KStest <- function(
  counts,
  classes,
  cores = 4
){
  
  #input checks
  if(any(grepl("-", classes))) {
    stop("The - character may not be included in the classes arg")
  }
  if(any(is.na(counts))) {
    stop("NAs in counts")
  }
  if(!is.numeric(counts)) {
    stop("non-numeric values in counts")
  }
  if(length(classes) != ncol(counts)) {
    stop("length(classes) != ncol(counts)")
  }
  if(nrow(counts) == 1) {
    stop("Only 1 gene in counts.")
  }
  
  #find unique classes to compare
  uc <- unique(classes)
  names(uc) <- uc
  
  #find all combinations of unique classes
  n <- combn(uc, 2)
  n <- apply(n, 2, sort)
  cmbNames <- paste(n[1, ], n[2, ], sep = "-")
  colnames(n) <- cmbNames
  n <- as.data.frame(n)
  
  #setup x variable for ks.test
  c1 <- splitCountsByClass(n, counts, classes, 1)
  
  #setup y variable for ks.test
  c2 <- splitCountsByClass(n, counts, classes, 2)
 
 #run ks.test for each class comparison (outer; mclapply) and
 # gene (inner; map2_dfr). Name combinations in object, reformat to a data.frame
 # and then tibble.
 parallel::mclapply(1:length(c1), function(i) {
   x <- c1[[i]]
   y <- c2[[i]]
   purrr::map2_dfr(x, y, ~runKS(.x, .y), .id = "gene")
 }, mc.cores = cores) %>%
  setNames(cmbNames) %>%
  dplyr::bind_rows(.id = "combination") %>%
  tibble::as_tibble()
}

#sets up x and y variables for ks.test
#n: A martix with each column containing a combination, of 2 elements, of the
# unique classes in the classes arg. Typically output from combn(classes, 2).
# Naming the colnames of the matrix facilitates since this function will return
# a named list.
#exp: See exp param specifications for the KStest function.
#classes: See classes param specifications for the KStest function.
#idx: the row index (either 1 or 2) in n. Indicates if either x or y arg to
# ks.test is being setup.

splitCountsByClass <- function(n, exp, classes, idx) {
  lapply(n, function(x) {
    as.data.frame(t(exp[, classes == x[idx]]))
  })
}

#accepts the x and y args to ks.test and returns the results in a data.frame
runKS <- function(z, w) {
  res <- my.ks.test(z, w)
  data.frame(
    stat = res[[1]],
    p.value = res[[2]]
  )
}

#' processKStest
#'
#' Identifies genes for classes that were significant in all comparisons using
#' the KStest function. In addition, calculates the sum of statistics for each
#' gene and cell type.
#'
#' @name processKStest
#' @rdname processKStest
#' @aliases processKStest
#' @param results Write something
#' @param classes Write something
#' @param alpha Write something
#' @return Write something
#' @author Jason T. Serviss
#' @examples
#' #nothing here yet
NULL

#' @export
#' @importFrom purrr map_dfr
#' @importFrom dplyr group_by summarize filter
#' @importFrom stringr str_replace

processKStest <- function(results, classes, alpha) {

  #find the unique classes
  uc <- unique(classes)
  names(uc) <- uc
  
  #process results
  purrr::map_dfr(uc, function(x, results) {
    .processResults(results, x, alpha)
  }, results = results, .id = "id")
}

#called in a loop over unique classes and accepts the mapped ks.test output and
# the current class (currClass). Summarizes the results by finding p.values that
# were significant in all comparisons. In addition, calculates the sum of the
# test statistics for each gene for the current class.
.processResults <- function(results, currClass, alpha) {
  results %>%
    dplyr::filter(
      stringr::str_replace(.data$combination, "(.*)-.*", "\\1") == currClass |
      stringr::str_replace(.data$combination, ".*-(.*)", "\\1") == currClass
    ) %>%
    dplyr::group_by(.data$gene) %>%
    dplyr::summarize(
      sigBool = all(.data$p.value < alpha),
      statSum = sum(.data$stat)
    )
}

