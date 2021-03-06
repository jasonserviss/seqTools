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

#' nTopDeltaCV
#'
#' Selected informative genes by fitting a support-vector regression to the
#' coefficient of variation (CV) as a function of the mean, and selecting genes
#' having the greatest offset from the fitted curve; this would correspond to
#' genes with higher-than-expected variance. Adpoted from Linnarsson lab python
#' implementation:
#' https://github.com/linnarsson-lab/cytograph/blob/master/cytograph/feature_selection.py
#'
#' @name nTopDeltaCV
#' @rdname nTopDeltaCV
#' @author Jason T. Serviss
#' @param cpm matrix; Matrix containing cpm values.
#' @param n Number of genes to select.
#' @keywords nTopDeltaCV
#'
#'
#' @export
NULL

#' @importFrom e1071 svm
#' @importFrom matrixStats rowMeans2 rowSds

nTopDeltaCV <- function(counts, n, rm = 0) {
  valid <- matrixStats::rowSums2(counts) > rm
  mu <- matrixStats::rowMeans2(counts)
  sd <- matrixStats::rowSds(counts)
  ok <- mu > 0 & sd > 0
  cv <- sd[ok] / mu[ok]
  
  log2_m <- log2(mu[ok])
  log2_cv <- log2(cv)
  
  svr_gamma <- 1000 / length(mu[ok])
  modelsvm <- svm(log2_cv ~ log2_m, gamma = svr_gamma)
  score <- log2_cv - predict(modelsvm, log2_m)
  score <- score * valid[ok]
  names(score) <- rownames(counts)[ok]
  s <- sort(score, decreasing = TRUE)[1:n]
  match(names(s), rownames(counts))
}

#' recursiveFE
#'
#' Recursive Feature Elimination; First, the algorithm fits the model to all
#' predictors. Each predictor is ranked using it’s importance to the model. Let
#' S be a sequence of ordered numbers which are candidate values for the number
#' of predictors to retain (S1 > S2, …). At each iteration of feature selection,
#' the Si top ranked predictors are retained, the model is refit and performance
#' is assessed. The value of Si with the best performance is determined and the
#' top Si predictors are used to fit the final model.
#'
#' @name recursiveFE
#' @rdname recursiveFE
#' @author Jason T. Serviss
#' @param cpm Matrix; Matrix containing cpm values.
#' @param classes Character; classifications corresponding to the cpm columns.
#' @param nCV Integer; Number of cross-fold validations.
#' @param s Numeric; Vector indicating the number of features to select in
#' each testing round.
#' @keywords recursiveFE
#'
#'
#' @export
NULL

#' @importFrom caret rfeControl rfe predictors rfFuncs

recursiveFE <- function(cpm, classes, nCV = 10, s) {
  normalized <- scale(cpm, center = TRUE, scale = TRUE)
  control <- rfeControl(functions = caret::rfFuncs, method = "cv", number = nCV)
  results <- rfe(
    as.data.frame(t(normalized)),
    as.factor(classes),
    sizes = s,
    rfeControl = control
  )
  p <- predictors(results)
  which(rownames(cpm) %in% p)
}

#' findHighlyCorrelated
#'
#' Identifies genes with a level of correlation greater than a specified cutoff.
#' Returns indexes of all genes that are NOT correlated > the specified cutoff.
#'
#' @name findHighlyCorrelated
#' @rdname findHighlyCorrelated
#' @author Jason T. Serviss
#' @param cpm Matrix; Matrix containing cpm values.
#' @param cut Integer; Genes with a correlation > cut will be reported.
#' @keywords recursiveFE
#' @export
NULL

#' @importFrom caret findCorrelation

findHighlyCorrelated <- function(cpm, cut) {
  correlationMatrix <- cor(t(cpm))
  hc <- findCorrelation(correlationMatrix, cutoff = cut)
  n <- 1:nrow(cpm)
  n[!n %in% hc]
}

#' foldChangePerClass
#'
#' Calculates fold change for each gene with each class vs all other classes.
#'
#' @name foldChangePerClass
#' @rdname foldChangePerClass
#' @author Jason T. Serviss
#' @param counts The matrix holding expression values.
#' @param class A tibble with columns \emph{class} and \emph{sample} indicating
#'    the class and sample ID respectivley.
#' @keywords foldChangePerClass
#'
#'
#' @export
#' @importFrom tibble tibble
#' @importFrom dplyr filter
#' @importFrom matrixStats rowMeans2

NULL

foldChangePerClass <- function(counts, classes) {
  uGroups <- unique(classes$class)
  
  res <- sapply(1:length(uGroups), function(x) {
    samplesA <- filter(classes, class == uGroups[x])$sample
    samplesB <- filter(classes, class != uGroups[x])$sample
    a <- rowMeans2(counts[, colnames(counts) %in% samplesA])
    b <- rowMeans2(counts[, colnames(counts) %in% samplesB])
    a/b
  })
  colnames(res) <- uGroups
  rownames(res) <- rownames(counts)
  return(res)
}

#' pearsonsCor
#'
#' Calculates 1-Pearson's correlation between columns and returns results as a
#' dist object.
#'
#' @name pearsonsCor
#' @rdname pearsonsCor
#' @author Jason T. Serviss
#' @param counts The matrix holding expression values.
#' @param select Optional; row indexes to include.
#' @keywords pearsonsCor
#'
#'
#' @export
NULL

pearsonsCor <- function(cpm, select = NULL) {
  if(is.null(select)) select <- 1:nrow(cpm)
  as.dist(1 - cor(cpm[select, ], method = "p"))
}
