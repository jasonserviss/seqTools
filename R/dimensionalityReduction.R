#' run.tsne
#'
#'
#' Calculates the x and y coordinates of the mean of each classified group.
#'
#'
#' This method is typically only used in conjunction with plotting. It
#' calculates the 2 dimensional location of the mean of each classified group
#' in the supplied unsupervised dimensionality reduction (t-SNE) data
#' representation.
#'
#' @name run.tsne
#' @rdname run.tsne
#' @aliases run.tsne
#' @param my.dist A distance object typically produced with pearsons.dist.
#' @param dims Argument to Rtsne. Numeric indicating the output dimensions.
#' @param theta Argument to
#'    [Rtsne](https://cran.r-project.org/web/packages/Rtsne/index.html).
#' @param initial_dims Argument to
#'    [Rtsne](https://cran.r-project.org/web/packages/Rtsne/index.html).
#' @param max_iter Argument to
#'    [Rtsne](https://cran.r-project.org/web/packages/Rtsne/index.html).
#' @param perplexity Argument to
#'    [Rtsne](https://cran.r-project.org/web/packages/Rtsne/index.html).
#' @param seed The desired seed to set before running.
#' @param is_distance Argument to
#'    [Rtsne](https://cran.r-project.org/web/packages/Rtsne/index.html).
#' @param ... Additional arguments to pass on
#' @return A matrix containing the tSNE dimensionally reduced data.
#' @author Jason T. Serviss
#' @keywords run.tsne
#' @examples
#' my.dist <- pearsons.dist(testingCounts[, 2:12], 1:14)
#' tsne <- run.tsne(my.dist, max_iter = 10, perplexity = 2)
#'
NULL

#' @rdname run.tsne
#' @importFrom Rtsne Rtsne
#' @export

run.tsne <- function(
  my.dist, dims = 2, theta = 0, initial_dims = 50, max_iter = 2000,
  perplexity = 10, seed = 11, is_distance = TRUE, ...
){
  set.seed(seed)
  
  my.tsne <- Rtsne(
    my.dist, dims = dims, initial_dims = initial_dims, max_iter = max_iter,
    perplexity = perplexity, theta = theta, is_distance = is_distance
  )$Y
  
  rownames(my.tsne) <- attr(my.dist, "Labels")
  return(my.tsne)
}

#' run.umap
#'
#' pak::pak("jlmelville/uwot")
#'
#'
#' @name run.umap
#' @rdname run.umap
#' @aliases run.umap
#' @param my.dist A distance object typically produced with pearsons.dist.
#' @param ... Additional arguments to pass on
#' @return
#' @author Jason T. Serviss
#' @keywords run.umap
#' @examples
#'
#' tsne <- run.umap(my.dist, max_iter = 10)
#'
NULL

#' @rdname run.umap
#' @importFrom uwot umap2
#' @export

run.umap <- function(
  cpm, n_neighbors = 30, min_dist = 0.3, metric = "euclidean", ...
){
  set.seed(seed)
  #alternative method to run umap on principal components?
  
  my.umap <- umap2(
    cpm.log2(cpm), n_neighbors = n_neighbors, min_dist = min_dist,
    metric = metric
  )
  
  rownames(my.umap) <- rownames(cpm)
  return(my.umap)
}

#UMAP
#uwot package, umap2 function

#PCA



#PCA

