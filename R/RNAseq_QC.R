#' detectLowQualityGenes
#'
#' In gene expression counts data if can often be the case that some genes are
#' not detected. This can simply be due to the fact that the gene is not
#' expressed in the tissue or, in addition, that the sequencing depth was not
#' sufficient to detect the gene. In addition, some genes may be detected but in
#' so few samples or at such a low level that it makes the quantified value
#' highly unreliable. In these cases, it is desireable to remove the gene before
#' downstream analysis which is facilitated by this function.
#'
#' @name detectLowQualityGenes
#' @rdname detectLowQualityGenes
#' @param counts data.frame; A data frame with counts data with gene names as
#' rownames.
#' @param mincount numeric; A minimum rowSum for which rows with a higher rowSum
#' will be detected. Default = 0.
#' @return A logical vector with length = nrow(counts) that is TRUE when the
#' counts data.frame row meets both tested conditions.
#' @author Jason Serviss
#' @examples
#'
#' counts <- data.frame(c(0, runif(10)), c(0, runif(10)), c(0, runif(10)))
#' detectLowQualityGenes(counts)
#'
NULL
#' @export

detectLowQualityGenes <- function(
  counts,
  mincount = 0
){
  #input checks

  bool <- rowSums(counts) > mincount

  message <- paste0(
    "Detected ", sum(!bool), " low quality genes out of ", nrow(counts),
    " genes input (", round(100 * (sum(!bool) / nrow(counts)), digits = 2),
    "%)."
  )
  print(message)
  return(bool)
}

#' detectLowQualityCells.totalCounts
#'
#' It is often the case that some samples from sequencing experiments are of
#' low quality, in many cases due to issues during the sample preperation stage.
#' Due to the fact that these samples represent a high level of technical noise,
#' it is often desirable to remove these before downstream analysis which is
#' facilitated by this function. The function achieves this by detecting cells
#' whose sum across all genes is > mincount.
#'
#' @name detectLowQualityCells.totalCounts
#' @rdname detectLowQualityCells.totalCounts
#' @param counts data.frame; A data frame with counts data with gene names as
#' rownames and sample names as colnames.
#' @param mincount numeric; A minimum colSum for which columns with a higher
#' colSum will be detected. Default = 4e5.
#' @return A logical vector with length = ncol(counts) that is TRUE when the
#' counts data.frame column contains a sample with colSums > mincount.
#' @author Jason Serviss
#' @examples
#' set.seed(8292)
#' x <- runif(2e4)
#' y <- runif(2e4, 1, 100)
#' names <- paste0(letters, 1:2e4)
#' counts <- data.frame(a = x, b = y, c = y, row.names = names)
#' detectLowQualityCells.totalCounts(counts)
#'
NULL
#' @export
#' @importFrom stats quantile

detectLowQualityCells.totalCounts <- function(
  counts,
  mincount = 4e5
){
  #setup output vector
  output <- vector(mode = "logical", length = ncol(counts))
  names(output) <- colnames(counts)

  #colsums check
  cs <- colSums(counts) > mincount
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

#' detectLowQualityCells.housekeeping
#'
#' It is often the case that some samples from sequencing experiments are of
#' low quality, in many cases due to issues during the sample preperation stage.
#' Due to the fact that these samples represent a high level of technical noise,
#' it is often desirable to remove these before downstream analysis which is
#' facilitated by this function. The function achieves this by utilizing a house
#' keeping gene and assuming its log2 expression to be normally distributed.
#' We then detect samples where the probability of the expression for the house
#' keeping gene in that sample is greater than the quantile.cut argument.
#'
#' @name detectLowQualityCells.housekeeping
#' @rdname detectLowQualityCells.housekeeping
#' @param counts data.frame; A data frame with counts data with gene names as
#' rownames and sample names as colnames.
#' @param geneName character; The gene name to use for the quantile cutoff. This
#' must be present in the rownames of the counts argument. Default is ACTB.
#' @param quantileCut numeric; This indicates probability at which the quantile
#' cutoff will be calculated using the normal distribution. Default = 0.01.
#' @return A logical vector with length = ncol(counts) that is TRUE when the
#' counts data.frame column contains a sample with meeting the criteria specified
#' by the arguments.
#' @author Jason Serviss
#' @examples
#' set.seed(8292)
#' x <- runif(2e4)
#' y <- runif(2e4, 1.5, 100)
#' names <- paste0(letters, 1:2e4)
#' counts <- data.frame(a = x, b = y, c = y, row.names = names)
#' detectLowQualityCells.housekeeping(counts, geneName = "a1")
#'
NULL
#' @export
#' @importFrom stats median qnorm

detectLowQualityCells.housekeeping <- function(
  counts,
  geneName = 'ACTB',
  quantileCut = 0.01
){
  #input checks
  ##check that geneName is in rownames counts
  if(!geneName %in% rownames(counts)) {
    stop("geneName is not found in rownames(counts)")
  }

  #setup output vector
  output <- vector(mode = "logical", length = ncol(counts))
  names(output) <- colnames(counts)

  #house keeping check
  counts.log <- .norm.log.counts(counts)
  cl.act <- counts.log[geneName, colSums(counts) != 0]
  cl.act.m <- median(cl.act)
  cl.act.sd <- sqrt(
    sum((cl.act[cl.act > cl.act.m] - cl.act.m) ^ 2) /
      (sum(cl.act > cl.act.m) - 1)
  )
  my.cut <- qnorm(p = quantileCut, mean = cl.act.m, sd = cl.act.sd)
  bool <- counts.log[geneName, ] >= my.cut
  output[bool] <- TRUE

  message <- paste0(
    "Detected ", sum(!output), " low quality cells out of ", ncol(counts),
    " cells input (", round(100 * (sum(!output) / ncol(counts)), digits = 2),
    "%) based on ", geneName, " expression."
  )
  print(message)
  return(output)
}

#calculates log cpm
.norm.log.counts <- function(counts) {
  log2(.norm.counts(counts) + 1)
}

#calculates cpm
.norm.counts <- function(counts) {
  t(t(counts) / colSums(counts) * 10^6)
}

#' detectLowQualityCells.ERCCfrac
#'
#' It is often the case that some samples from sequencing experiments are of
#' low quality, in many cases due to issues during the sample preperation stage.
#' Due to the fact that these samples represent a high level of technical noise,
#' it is often desirable to remove these before downstream analysis which is
#' facilitated by this function. The function achieves this by indicating samples
#' that have 0 ERCC reads and those that have a fraction of ERCC reads .
#'
#' @name detectLowQualityCells.ERCCfrac
#' @rdname detectLowQualityCells.ERCCfrac
#' @param counts data.frame; A data frame with counts data with gene names as
#' rownames and sample names as colnames.
#' @param geneName character; The gene name to use for the quantile cutoff. This
#' must be present in the rownames of the counts argument. Default is ACTB.
#' @param quantileCut numeric; This indicates probability at which the quantile
#' cutoff will be calculated using the normal distribution. Default = 0.99.
#' @return A logical vector with length = ncol(counts) that is TRUE when the
#' counts data.frame column contains a sample with meeting the criteria specified
#' by the arguments.
#' @author Jason Serviss
#' @examples
#' set.seed(8292)
#' x <- runif(2e4)
#' y <- runif(2e4, 1.5, 100)
#' names <- paste0(letters, 1:2e4)
#' counts <- data.frame(a = x, b = y, c = y, row.names = names)
#' detectLowQualityCells.housekeeping(counts, geneName = "a1")
#'
NULL
#' @export

detectLowQualityCells.ERCCfrac <- function(
  counts,
  ercc,
  quantileCut = 0.99
){
  #setup output vector
  output <- vector(mode = "logical", length = ncol(counts))
  names(output) <- colnames(counts)

  #calculate fraction of ercc
  cs <- colSums(counts)
  cs.ercc <- colSums(ercc)
  frac.ercc <-  cs.ercc / (cs.ercc + cs)
  frac.ercc[is.na(frac.ercc)] <- 0 #when cs == 0
  l.frac.ercc <- log2(frac.ercc)

  #calculate normal approximation and cut
  #median less sensitive to outliers
  #is.infinite needed for when cs.ercc == 0
  m <- median(l.frac.ercc[!is.infinite(l.frac.ercc)], na.rm = TRUE)
  sd <- sd(l.frac.ercc[!is.infinite(l.frac.ercc)], na.rm = TRUE)
  p.cut.h <- 2^qnorm(quantileCut, m, sd)
  p.cut.l <- 2^qnorm(1 - quantileCut, m, sd)
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
