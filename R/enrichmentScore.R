
#' enrichmentScore
#'
#' Calculates an the enrichment of genes per classification.
#' https://www.biorxiv.org/content/biorxiv/early/2018/04/06/294918.full.pdf
#'
#'
#' @name enrichmentScore
#' @rdname enrichmentScore
#' @aliases enrichmentScore
#' @param cpm matrix; Matrix containing cpm values.
#' @param classes Character; Classes corresponding to columns in cpm argument.
#' @return Returns the a list of enrichment scores with one element per 
#' classification.
#' @author Jason T. Serviss
#' @keywords nTopVar
#'
NULL

#' @rdname enrichmentScore
#' @importFrom matrixStats rowMeans2
#' @importFrom future.apply future_lapply
#' @export

enrichmentScore <- function(cpm, classes, e1 = 0.1, e2 = 0.01) {
  
  if(!is.character(classes)) stop("classes arg should be a character vector")
  options(future.globals.maxSize = +Inf)
  
  uClasses <- unique(classes)
  scores <- future_lapply(1:length(uClasses), FUN = function(i) {
    currClass <- uClasses[i]
    
    #get cluster samples
    cGenes <- cpm[, classes == currClass]
    
    #get non cluster samples
    ncGenes <- cpm[, classes != currClass]
    
    #calculate fractions
    f <- apply(cGenes, 1, function(x) length(which(x != 0)))
    fj <- apply(ncGenes, 1, function(x) length(which(x != 0)))
    
    #calculate means
    u <- matrixStats::rowMeans2(cGenes)
    uj <- matrixStats::rowMeans2(ncGenes)
    
    #calculate enrichment
    e <- ((f + e1) / (fj + e1)) * ((u + e2) / (uj + e2))
    e <- sort(e, decreasing = TRUE)
    e
  })
  names(scores) <- uClasses
  scores
}
