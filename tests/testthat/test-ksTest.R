context("test-kstest.R")

################################################################################
#                                                                              #
# my.ks.test                                                                   #
#                                                                              #
################################################################################

test_that("my.ks.test without ties, args equal length", {
  
  #prepare input data without ties
  set.seed(29200)
  test <- TRUE
  while(test) {
    x <- rnorm(10000)
    y <- rnorm(10000)
    test <- any(identical(x, y))
  }
  
  #setup expected
  res <- ks.test(x, y, alternative = "two.sided", exact = TRUE)
  exp_s <- unname(res[[1]])
  exp_p <- unname(res[[2]])
  
  #run function
  output <- my.ks.test(x, y)
  out_s <- output[[1]]
  out_p <- output[[2]]
  
  #test
  expect_identical(exp_s, out_s)
  expect_identical(exp_p, out_p)
})

test_that("my.ks.test with ties, args equal length", {
  
  #prepare input data with ties
  x <- 1:100
  y <- 51:150
  
  #setup expected
  #ks.test will warn about ties which causes testthat to detect a warning
  res <- expect_warning(ks.test(x, y))
  exp_s <- unname(res[[1]])
  exp_p <- unname(res[[2]])
  
  #run function
  output <- my.ks.test(x, y, alternative = "two.sided", exact = TRUE)
  out_s <- output[[1]]
  out_p <- output[[2]]
  
  #test
  expect_identical(exp_s, out_s)
  expect_identical(exp_p, out_p)
})

test_that("my.ks.test without ties, args unequal length", {
  
  #prepare input data without ties
  set.seed(29200)
  test <- TRUE
  while(test) {
    x <- rnorm(100)
    y <- rnorm(50)
    test <- any(identical(x, y))
  }
  
  #setup expected
  res <- ks.test(x, y, alternative = "two.sided", exact = TRUE)
  exp_s <- unname(res[[1]])
  exp_p <- unname(res[[2]])
  
  #run function
  output <- my.ks.test(x, y)
  out_s <- output[[1]]
  out_p <- output[[2]]
  
  #test
  expect_identical(exp_s, out_s)
  expect_identical(exp_p, out_p)
})

test_that("my.ks.test with ties, args unequal length", {
  
  #prepare input data with ties
  x <- 1:100
  y <- 1:50
  
  #setup expected
  #ks.test will warn about ties which causes testthat to detect a warning
  expect_warning(res <- ks.test(x, y, alternative = "two.sided", exact = TRUE))
  exp_s <- unname(res[[1]])
  exp_p <- unname(res[[2]])
  
  #run function
  output <- my.ks.test(x, y)
  out_s <- output[[1]]
  out_p <- output[[2]]
  
  #test
  expect_identical(exp_s, out_s)
  expect_identical(exp_p, out_p)
})

test_that("my.ks.test input check NAs in x", {
  
  #prepare input
  x <- c(1:10, NA)
  y <- 1:10
  
  #test
  expect_error(my.ks.test(x, y, check = TRUE))
})

test_that("my.ks.test input check NAs in y", {
  
  #prepare input
  x <- 1:10
  y <- c(1:10, NA)
  
  #test
  expect_error(my.ks.test(x, y, check = TRUE))
})

test_that("my.ks.test input check !is.numeric(x)", {
  
  #prepare input
  x <- letters[1:10]
  y <- 1:10
  
  #test
  expect_error(my.ks.test(x, y, check = TRUE))
})

test_that("my.ks.test input check !is.numeric(y)", {
  
  #prepare input
  x <- 1:10
  y <- letters[1:10]
  
  #test
  expect_error(my.ks.test(x, y, check = TRUE))
})

test_that("my.ks.test input check length(x) > 1L", {
  
  #prepare input
  x <- numeric()
  y <- 1:10
  
  #test
  expect_error(my.ks.test(x, y, check = TRUE))
})

test_that("my.ks.test input check length(y) > 1L", {
  
  #prepare input
  x <- 1:10
  y <- numeric()
  
  #test
  expect_error(my.ks.test(x, y, check = TRUE))
})

################################################################################
#                                                                              #
# splitCountsByClass                                                           #
#                                                                              #
################################################################################

test_that("splitCountsByClass returns correct result with multiple genes", {
  
  #prepare input data
  exp <- matrix(rep(1:4, each = 10), nrow = 2, dimnames = list(letters[1:2], NULL))
  classes <- rep(LETTERS[1:4], each = 5)
  idx <- 1

  #setup n arg
  uc <- unique(classes)
  names(uc) <- uc
  n <- combn(uc, 2)
  cmbNames <- paste(n[1, ], n[2, ], sep = "-")
  colnames(n) <- cmbNames
  n <- as.data.frame(n)
  
  #setup expected
  A <- data.frame(a = rep(1L, 5), b = rep(1L, 5))
  B <- data.frame(a = rep(2L, 5), b = rep(2L, 5))
  C <- data.frame(a = rep(3L, 5), b = rep(3L, 5))
  expected <- list(A, A, A, B, B, C)
  names(expected) <- cmbNames
  
  #run function
  output <- splitCountsByClass(n, exp, classes, idx)
  
  #test
  expect_identical(expected, output)
})

################################################################################
#                                                                              #
# runKS                                                                        #
#                                                                              #
################################################################################

test_that("runKS returns correct result", {
  
  #this should be covered by testing of the my.ks.test function
  #prepare input data
  set.seed(29200)
  test <- TRUE
  while(test) {
    x <- rnorm(10)
    y <- rnorm(10)
    test <- any(identical(x, y))
  }
  #setup expected
  test <- ks.test(x, y)
  e.stat <- unname(test[[1]])
  e.p <- unname(test[[2]])

  #run function
  output <- runKS(x, y)
  out.stat <- output[[1]]
  out.p <- output[[2]]
  
  #test
  expect_identical(e.stat, out.stat)
  expect_identical(e.p, out.p)
})

################################################################################
#                                                                              #
# KStest                                                                       #
#                                                                              #
################################################################################

test_that("KStest errors with - character present in classes arg", {
  
  #prepare input data
  exp <- matrix(1:10, ncol = 2)
  classes <- c("A-B", "B-C")
  
  #test
  expect_error(KStest(exp, classes, cores = 1))
})

test_that("KStest errors with non-numeric in counts", {
  
  #prepare input data
  exp <- matrix(c(letters, 1:10), ncol = 2)
  classes <- LETTERS[1:2]
  
  #test
  expect_error(KStest(exp, classes, cores = 1))
})

test_that("KStest errors with NA in counts", {
  
  #prepare input data
  exp <- matrix(c(NA, 1:9), ncol = 2)
  classes <- LETTERS[1:2]

  #test
  expect_error(KStest(exp, classes, cores = 1))
})

test_that("KStest errors without matching dims", {
  
  #prepare input data
  exp <- matrix(c(NA, 1:9), ncol = 2)
  classes <- LETTERS[1:5]
  
  #test
  expect_error(KStest(exp, classes, cores = 1))
})

test_that("KStest errors nrow(counts) = 1", {
  
  #prepare input data
  exp <- matrix(c(NA, 1:9), nrow = 1)
  classes <- LETTERS[1:5]
  
  #test
  expect_error(KStest(exp, classes, cores = 1))
})

test_that("KStest returns the correct result format", {
  
  #prepare input data
  rows <- 100
  cols <- 10
  exp <- matrix(
    1:1000,
    ncol = cols,
    nrow = rows,
    dimnames = list(
    paste(letters, 1:rows, sep = ""),
      paste(LETTERS[1:cols], 1:cols, sep = ""))
  )
  classes <- rep(LETTERS[1:2], each = 5)
  cores <- 1
  
  #setup expected
  nr <- nrow(exp) * ncol(combn(unique(classes), 2))
  nc <- 4L
  cnames <- c("combination", "gene", "stat", "p.value")
  c <- c("tbl_df", "tbl", "data.frame")
  c1 <- "character"
  c2 <- "character"
  c3 <- "numeric"
  c4 <- "numeric"
  uGenes <- rownames(exp)
  nCombs <- ncol(combn(unique(classes), 2))
  combs <- paste(
    combn(unique(classes), 2)[1, ],
    combn(unique(classes), 2)[2, ],
    sep = "-"
  )
  combTotal <- rep(combs, each = nrow(exp))
  genesTotal <- rep(rownames(exp), nCombs)
  
  #run function
  output <- KStest(exp, classes, cores)
  
  #test
  expect_identical(nr, nrow(output)) #correct dimensions
  expect_identical(nc, ncol(output)) #correct dimensions
  expect_identical(cnames, colnames(output)) #correct colnames
  expect_identical(c, class(output)) #correct class
  expect_identical(c1, class(output[[1]])) #correct class col 1
  expect_identical(c2, class(output[[2]])) #correct class col 2
  expect_identical(c3, class(output[[3]])) #correct class col 3
  expect_identical(c4, class(output[[4]])) #correct class col 4
  expect_identical(unique(output[[2]]), uGenes) #all unique genes present
  expect_identical(unique(output[[1]]), combs) #all unique combinations present
  expect_identical(genesTotal, output[[2]]) #genes freq
  expect_identical(combTotal, output[[1]]) #combination freq
})

test_that("KStest returns the correct result content", {
  
  #prepare input data
  set.seed(78340923)
  x1 <- rnorm(10, 1, 1)
  y1 <- rpois(10, 1)
  
  exp <- matrix(
    c(x1, x1, x1, y1),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(letters[1:2], NULL)
  )
  classes <- rep(LETTERS[1:2], each = 10)
  cores <- 1
  
  #setup expected
  expected <- tibble::tibble(
    combination = c("A-B", "A-B"),
    gene = letters[1:2],
    stat = expect_warning(c(ks.test(x1, x1)[[1]], ks.test(x1, y1)[[1]])),
    p.value = expect_warning(c(ks.test(x1, x1)[[2]], ks.test(x1, y1)[[2]]))
  )
  
  #run function
  output <- KStest(exp, classes, cores)
  
  #test
  expect_identical(expected, output)
})

################################################################################
#                                                                              #
# processKStest                                                                #
#                                                                              #
################################################################################

test_that("processKStest returns the correct result", {
  
  #prepare input data
  set.seed(78340923)
  x1 <- rnorm(100, 1, 1)
  y1 <- rpois(100, 1)
  
  exp <- matrix(
    c(x1, y1, y1, y1, rep(x1, 4)),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(letters[1:2], NULL)
  )
  classes <- rep(LETTERS[1:4], each = 100)
  cores <- 1
  input <- KStest(exp, classes, cores)
  
  #setup expected
  statSum <- .27 * 3
  
  #run function
  output <- processKStest(input, classes, 0.05)
  
  #test
  expect_true(filter(output, id == "A" & gene == "a")[[3]])
  expect_true(all(filter(output, id != "A" | gene == "b")[[3]] == FALSE))
  expect_identical(filter(output, id == "A" & gene == "a")[[4]], statSum)
  expect_true(all(filter(output, id != "A" | gene == "b")[[4]] < statSum))
})
