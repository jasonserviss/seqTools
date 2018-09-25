JackRandom <- function(
  scaled.data,
  prop.use = 0.01,
  r1.use = 1,
  r2.use = 5,
  seed.use = 1,
  rev.pca = FALSE,
  weight.by.var = weight.by.var,
  maxit = 1000
) {
  set.seed(seed = seed.use)
  rand.genes <- sample(
    x = rownames(x = scaled.data),
    size = nrow(x = scaled.data) * prop.use
  )
  # make sure that rand.genes is at least 3
  if (length(x = rand.genes) < 3){
    rand.genes <- sample(x = rownames(x = scaled.data), size = 3)
  }
  data.mod <- scaled.data
  data.mod[rand.genes, ] <- MatrixRowShuffle(x = scaled.data[rand.genes, ])
  temp.object <- new("seurat")
  temp.object@cell.names <- colnames(x = data.mod)
  temp.object@scale.data <- data.mod
  temp.object <- RunPCA(
    object = temp.object,
    pcs.compute = r2.use,
    pc.genes = rownames(x = data.mod),
    rev.pca = rev.pca,
    weight.by.var = weight.by.var,
    do.print = FALSE,
    maxit = maxit
  )
  fake.x <- PCALoad(object = temp.object)
  fake.rot <- PCAEmbed(object = temp.object)
  return(fake.x[rand.genes, r1.use:r2.use])
}

#' Independently shuffle values within each row of a matrix
#'
#' Creates a matrix where correlation structure has been removed, but overall values are the same
#'
#' @param x Matrix to shuffle
#'
#' @return Returns a scrambled matrix, where each row is shuffled independently
#'
#' @importFrom stats runif
#'
#' @export
#'
#' @examples
#' mat <- matrix(data = rbinom(n = 25, size = 20, prob = 0.2 ), nrow = 5)
#' mat
#' MatrixRowShuffle(x = mat)
#'
MatrixRowShuffle <- function(x) {
  x2 <- x
  x2 <- t(x = x)
  ind <- order(c(col(x = x2)), runif(n = length(x = x2)))
  x2 <- matrix(
    data = x2[ind],
    nrow = nrow(x = x),
    ncol = ncol(x = x),
    byrow = TRUE
  )
  return(x2)
}