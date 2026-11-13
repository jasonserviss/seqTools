#' filter.counts.data
#'
#' This is a wrapper function for \code{\link{detect.low.quality.genes}},
#' \code{\link{detect.lqc_total.counts}},
#' \code{\link{detect.lqc_ERCC.frac}},
#' \code{\link{detect.lqc_housekeeping}}
#'
#' @name filter.counts.data
#' @rdname filter.counts.data
#' @param counts data.frame; The counts data.
#' @param counts.ERCC data.frame; The counts ERCC data.
#' @param filters Character; A vector indicating the types of filtering to be
#' performed. Options include: "genes", "total.counts",
#' "ERCC.frac", and "housekeeping".
#' @param gene.min.count See \code{\link{detect.low.quality.genes}}
#' @param cell.min.count See \code{\link{detect.lqc_total.counts}}
#' @param gene.name See \code{\link{detect.lqc_housekeeping}}
#' @param quantile.cut.hk See \code{\link{detect.lqc_housekeeping}}
#' @param quantile.cut.ercc See \code{\link{detect.lqc_ERCC.frac}}
#' @param percentile See \code{\link{detect.lqc_ERCC.frac}}.
#' @return A list with the filtered counts as the first element and the filtered
#' ERCC genes as the second element.
#' @author Jason Serviss
#'
NULL
#' @export

filter.counts.data <- function(
    counts, counts.ERCC,
    filters = c("genes", "total.counts", "ERCC.frac", "housekeeping"),
    gene.min.count, cell.min.count, gene.name, quantile.cut.hk,
    quantile.cut.ercc
){
  #remove low quality genes
  if("genes" %in% filters) {
    c <- counts[detect.low.quality.genes(counts, min.count = gene.min.count), ]
  }

  #remove low quality cells
  lqc.total.c <- lqc.housekeeping <- lqc.ERCC.frac <- rep(TRUE, ncol(counts))

  if("total.counts" %in% filters) {
    lqc.total.c <- detect.lqc_total.counts(counts, min.count = cell.min.count)
  }

  if("housekeeping" %in% filters) {
    lqc.housekeeping <- detect.lqc_housekeeping(counts,
      gene.name = gene.name, quantile.cut = quantile.cut.hk
    )
  }

  if("ERCCfrac" %in% filters) {
    lqc.ERCC.frac <- detect.lqc_ERCC.frac(counts, countsERCC,
      quantile.cut = quantile.cut.ercc
    )
  }

  lqc <- lqc.total.c & lqc.housekeeping & lqc.ERCC.frac

  print(paste0(
    "Removing a total of ",
    sum(!lqc),
    " cells based on the calculated metrics."
  ))

  c <- c[, lqc]
  counts.ERCC <- counts.ERCC[, lqc]

  return(list(c, counts.ERCC))
}


#' detect.low.quality.genes
#'
#' In gene expression counts data, it can often be the case that some genes are
#' not detected. This can simply be due to the fact that the gene is not
#' expressed in the tissue or that the sequencing depth was not
#' sufficient to detect the gene. In addition, some genes may be detected but in
#' so few samples or at such a low level that it makes the quantified value
#' highly unreliable. In these cases, it is desirable to remove the gene before
#' downstream analysis which is facilitated by this function.
#'
#' @name detect.low.quality.genes
#' @rdname detect.low.quality.genes
#' @param counts data.frame; A data frame with counts data with gene names as
#' rownames.
#' @param mincount numeric; A minimum rowSum for which rows with a higher rowSum
#' will be detected. Default = 0.
#' @return A logical vector with length = nrow(counts) that is TRUE when the
#' counts data.frame row meets the tested conditions.
#' @author Jason Serviss
#' @examples
#'
#' counts <- data.frame(c(0, runif(10)), c(0, runif(10)), c(0, runif(10)))
#' detect.low.quality.genes(counts)
#'
NULL
#' @export

detect.low.quality.genes <- function(
  counts,
  min.count = 0
){
  #input checks
  if(!is.data.frame(counts)) {stop("The counts vairable is not a data.frame.")}
  if(!is.numeric(min.count)) {stop("The mincount variable is not numeric.")}

  bool <- rowSums(counts) > min.count

  message <- paste0(
    "Detected ", sum(!bool), " low quality genes out of ", nrow(counts),
    " genes input (", round(100 * (sum(!bool) / nrow(counts)), digits = 2),
    "%)."
  )
  print(message)
  return(bool)
}

#' detect.lqc_total.counts
#'
#' It is often the case that some samples from sequencing experiments are of
#' low quality cells (lqc), in many cases due to issues during the sample
#' preparation stage.
#' Due to the fact that these samples represent a high level of technical noise,
#' it is often desirable to remove these before downstream analysis which is
#' facilitated by this function. The function achieves this by detecting cells
#' whose sum across all genes is > mincount.
#'
#' @name detect.lqc_total.counts
#' @rdname detect.lqc_total.counts
#' @param counts data.frame; A data frame with counts data with gene names as
#' rownames and sample names as colnames.
#' @param min.count numeric; A minimum colSum for which columns with a higher
#' colSum will be detected. Default = 4e5.
#' @return A logical vector with length = ncol(counts) that is TRUE where the
#' counts data.frame column contains a sample with colSums > mincount.
#' @author Jason Serviss
#' @examples
#' set.seed(8292)
#' x <- runif(2e4)
#' y <- runif(2e4, 1, 100)
#' names <- paste0(letters, 1:2e4)
#' counts <- data.frame(a = x, b = y, c = y, row.names = names)
#' detect.lqc_total.counts(counts)
#'
NULL
#' @export

detect.lqc_total.counts <- function(
  counts,
  min.count = 4e5
){
  #input checks
  if(!is.data.frame(counts)) {stop("The counts vairable is not a data.frame.")}
  if(!is.numeric(min.count)) {stop("The mincount variable is not numeric.")}

  #setup output vector
  output <- vector(mode = "logical", length = ncol(counts))
  names(output) <- colnames(counts)

  #colsums check
  cs <- colSums(counts) > min.count
  output[cs] <- TRUE

  if(sum(cs) < 2) {
    stop("One or less samples passed the colSums check.")
  }

  message <- paste0(
    "Detected ", sum(!output), " low quality cells out of ", ncol(counts),
    " cells input (", round(100 * (sum(!output) / ncol(counts)), digits = 2),
    "%) based on total counts."
  )
  print(message)
  return(output)
}

#' detect.lqc_housekeeping
#'
#' It is often the case that some samples from sequencing experiments are of
#' low quality cells (lqc), in many cases due to issues during the sample
#' preparation stage.
#' Due to the fact that these samples represent a high level of technical noise,
#' it is often desirable to remove these before downstream analysis which is
#' facilitated by this function. The function achieves this by utilizing a house
#' keeping gene and assuming its log2 expression to be normally distributed.
#' We then detect samples where the probability of the expression for the house
#' keeping gene in that sample is greater than the quantile.cut argument.
#'
#' @name detect.lqc_housekeeping
#' @rdname detect.lqc_housekeeping
#' @param counts data.frame; A data frame with counts data with gene names as
#' rownames and sample names as colnames.
#' @param gene.name character; The gene name to use for the quantile cutoff.
#' This must be present in the rownames of the counts argument. Default is ACTB.
#' @param quantile.cut numeric; This indicates probability at which the quantile
#' cutoff will be calculated using the normal distribution. Default = 0.01.
#' @return A logical vector with length = ncol(counts) that is TRUE when the
#' counts data.frame column contains a sample with meeting the criteria
#' specified by the arguments.
#' @author Jason Serviss
#' @examples
#' set.seed(8292)
#' x <- runif(2e4)
#' y <- runif(2e4, 1.5, 100)
#' names <- paste0(letters, 1:2e4)
#' counts <- data.frame(a = x, b = y, c = y, row.names = names)
#' detect.lqc_housekeeping(counts, geneName = "a1")
#'
NULL
#' @export
#' @importFrom stats median qnorm

detect.lqc_housekeeping <- function(
  counts,
  gene.name = 'ACTB',
  quantile.cut = 0.01
){
  #input checks
  if(!is.data.frame(counts)) {stop("The counts vairable is not a data.frame.")}
  if(!is.numeric(quantile.cut)) {
    stop("The quantileCut variable is not numeric.")
  }
  if(!gene.name %in% rownames(counts)) {
    stop("gene.name is not found in rownames(counts)")
  }

  #setup output vector
  output <- vector(mode = "logical", length = ncol(counts))
  names(output) <- colnames(counts)

  #house keeping check
  counts.log <- cpm.log2(counts)
  cl.act <- counts.log[gene.name, colSums(counts) != 0]

  #calculate median and sd
  cl.act.m <- median(cl.act)
  numerator <- sum((cl.act[cl.act > cl.act.m] - cl.act.m) ^ 2)
  denominator <- sum(cl.act > cl.act.m) - 1
  cl.act.sd <- sqrt(numerator / denominator)

  #establish cutoff
  my.cut <- qnorm(p = quantile.cut, mean = cl.act.m, sd = cl.act.sd)
  bool <- counts.log[gene.name, ] >= my.cut
  output[bool] <- TRUE

  message <- paste0(
    "Detected ", sum(!output), " low quality cells out of ", ncol(counts),
    " cells input (", round(100 * (sum(!output) / ncol(counts)), digits = 2),
    "%) based on ", gene.name, " expression."
  )
  print(message)
  return(output)
}

#' detect.lqc_ERCC.frac
#'
#' It is often the case that some samples from sequencing experiments are of
#' low quality cells (lqc), in many cases due to issues during the sample
#' preparation stage.
#' Due to the fact that these samples represent a high level of technical noise,
#' it is often desirable to remove these before downstream analysis which is
#' facilitated by this function. The function achieves this by indicating
#' samples that have 0 ERCC reads and those that have a low fraction of ERCC
#' reads.
#'
#' @name detect.lqc_ERCC.frac
#' @rdname detect.lqc_ERCC.frac
#' @param counts data.frame; A data frame with counts data with gene names as
#' rownames and sample names as colnames.
#' @param ercc data.frame; A data frame with the ercc counts data with ercc
#' names as rownames and sample names as colnames.
#' @param quantile.cut numeric; This indicates probability at which the quantile
#' cutoff will be calculated using the normal distribution. Default = 0.99.
#' @return A logical vector with length = ncol(counts) that is TRUE when the
#' counts data.frame column contains a sample with meeting the criteria
#' specified by the arguments.
#' @author Jason Serviss
#' @examples
#' set.seed(8292)
#' x <- runif(2e4)
#' y <- runif(2e4, 1.5, 100)
#' names <- paste0(letters, 1:2e4)
#' counts <- data.frame(a = x, b = y, c = y, row.names = names)
#' detect.lqc_ERCC.frac(counts, ercc = )
#'
NULL
#' @export

detect.lqc_ERCC.frac <- function(
  counts,
  ercc,
  quantile.cut = 0.99
){
  #INPUT CHECKS
  if(!is.data.frame(counts)) {stop("The counts vairable is not a data.frame.")}
  if(!is.numeric(quantile.cut)) {
    stop("The quantileCut variable is not numeric.")
  }

  #setup output vector
  output <- vector(mode = "logical", length = ncol(counts))
  names(output) <- colnames(counts)

  #calculate fraction of ercc
  cs <- colSums(counts)
  cs.ercc <- colSums(ercc)
  frac.ercc <-  cs.ercc / (cs.ercc + cs)
  frac.ercc[is.na(frac.ercc)] <- 0 #when cs == 0
  l.frac.ercc <- log2(frac.ercc)

  #calculate normal approximation and cutoff
  #median less sensitive to outliers
  #is.infinite needed for when cs.ercc == 0
  m <- median(l.frac.ercc[!is.infinite(l.frac.ercc)], na.rm = TRUE)
  sd <- sd(l.frac.ercc[!is.infinite(l.frac.ercc)], na.rm = TRUE)
  p.cut.h <- 2^qnorm(quantile.cut, m, sd)
  p.cut.l <- 2^qnorm(1 - quantile.cut, m, sd)
  bool.h <- frac.ercc <= p.cut.h
  bool.l <- frac.ercc >= p.cut.l
  output[bool.h & bool.l] <- TRUE

  message <- paste0(
    "Detected ", sum(!output), " low quality cells out of ", ncol(counts),
    " cells input (", round(100 * (sum(!output) / ncol(counts)), digits = 2),
    "%) based on ERCC fraction."
  )
  print(message)
  return(output)
}
