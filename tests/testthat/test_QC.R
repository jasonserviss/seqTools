
data(testingCounts)

test_that("check that the filter.counts.data function returns the correct output with expected data", {

  ##prepare input where the list elements and the vector elements are named
  counts <- testingCounts[1:13, ]
  rownames(counts) <- counts[[1]]
  counts[[1]] <- NULL
  
  #generate ercc counts
  set.seed(8923023)
  synth <- rnbinom(5 * 11, mu = 2^runif(5, 0, 5), size = 1)
  ercc <- matrix(as.numeric(synth), nrow = 5)
  
  #setup expected data
  expected <- list(
    counts[, -5],
    ercc[, -5]
  )
  
  ##run function
  output <- filter.counts.data(counts, ercc,
    gene.min.count = 0,
    cell.min.count = 1,
    gene.name = "ACTB",
    quantile.cut.hk = 0.01,
    quantile.cut.ercc = 0.99
  )

  ##test
  expect_identical(expected, output)
})

test_that("check that the detect.low.quality.genes function returns the correct output with expected data", {

  ##prepare input where the list elements and the vector elements are named
  counts <- testingCounts[1:13, ]
  rownames(counts) <- counts[[1]]
  counts[[1]] <- NULL
  
  mc <- 50
  
  #setup expected data
  n <- c("ACTB", letters[1:11], "ERCC-1")
  expected <- c(rep(TRUE, 3), FALSE, TRUE, FALSE, rep(TRUE, 5), rep(FALSE, 2))
  names(expected) <- n
  
  ##run function
  output <- detect.low.quality.genes(counts, mc)

  ##test
  expect_identical(expected, output)
})

test_that("check that the detect.lqc_total.counts function returns the correct output with expected data", {

  ##prepare input where the list elements and the vector elements are named
  counts <- testingCounts[1:13, ]
  rownames(counts) <- counts[[1]]
  counts[[1]] <- NULL
  colnames(counts) <- LETTERS[1:11]
  
  mc <- 130
  
  #setup expected data
  expected <- c(TRUE, FALSE, TRUE, FALSE, rep(TRUE, 6), FALSE)
  names(expected) <- LETTERS[1:11]
  
  ##run function
  output <- detect.lqc_total.counts(counts, mc)

  ##test
  expect_identical(expected, output)
})

test_that("check that the detect.lqc_housekeeping function returns the correct output with expected data", {

  ##prepare input where the list elements and the vector elements are named
  counts <- testingCounts[1:13, ]
  rownames(counts) <- counts[[1]]
  counts[[1]] <- NULL
  colnames(counts) <- LETTERS[1:11]
  
  gene.name <- 'ACTB'
  quantile.cut <- 0.01
  
  #setup expected data
  expected <- c(rep(TRUE, 4), FALSE, rep(TRUE, 6))
  names(expected) <- LETTERS[1:11]
  
  ##run function
  output <- detect.lqc_housekeeping(counts, gene.name, quantile.cut)

  ##test
  expect_identical(expected, output)
})

test_that("check that the detect.lqc_ERCC.frac function returns the correct output with expected data", {

  ##prepare input where the list elements and the vector elements are named
  counts <- testingCounts[1:13, ]
  rownames(counts) <- counts[[1]]
  counts[[1]] <- NULL
  colnames(counts) <- LETTERS[1:11]
  
  #generate ercc counts
  set.seed(8923023)
  synth <- rnbinom(5 * 11, mu = 2^runif(5, 0, 5), size = 1)
  ercc <- matrix(as.numeric(synth), nrow = 5)
  
  quantile.cut <- 0.99
  
  #setup expected data
  expected <- c(FALSE, rep(TRUE, 10))
  names(expected) <- LETTERS[1:11]
  
  ##run function
  output <- detect.lqc_ERCC.frac(counts, ercc, quantile.cut)

  ##test
  expect_identical(expected, output)
})
