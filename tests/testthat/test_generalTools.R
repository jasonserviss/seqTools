
test_that("check that the namedListToTibble function returns the correct output with expected data 1", {

  ##prepare input where the list elements and the vector elements are named
  a <- LETTERS[1:10]
  na <- rev(letters)[1:10]
  names(a) <- na

  b <- letters[1:10]
  nb <- paste(10:1, letters[1:10], sep = "")
  names(b) <- nb

  l <- list(a100=a, b100=b)

  #setup expected data
  expected <- tibble(
    names = c(rep("a100", 10), rep("b100", 10)),
    inner.names = c(na, nb),
    variables = c(LETTERS[1:10], letters[1:10])
  )

  ##run function
  output <- namedListToTibble(l)

  ##test
  expect_identical(expected, output)
})

test_that("check that the namedListToTibble function returns the correct output with expected data 2", {

  ##prepare input where only the list elements are named
  l <- list(a=LETTERS[1:10], b=letters[1:5])

  #setup expected data
  expected <- tibble(
    names = c(rep("a", 10), rep("b", 5)),
    variables = c(LETTERS[1:10], letters[1:5])
  )

  ##run function
  output <- namedListToTibble(l)

  ##test
  expect_identical(expected, output)
})

test_that("check that the namedListToTibble function returns the correct output NA in vector name", {

  ##prepare input where only the list elements are named
  a <- LETTERS[1:10]
  na <- rev(letters)[1:5]
  names(a) <- na

  b <- letters[1:10]
  nb <- paste(10:1, letters[1:10], sep = "")
  names(b) <- nb

  l <- list(a100=a, b100=b)

  #setup expected data
  expected <- tibble(
    names = c(rep("a100", 10), rep("b100", 10)),
    inner.names = c(na, rep(NA, 5), nb),
    variables = c(LETTERS[1:10], letters[1:10])
  )

  ##run function
  output <- namedListToTibble(l)

  ##test
  expect_identical(expected, output)
})

test_that("check that the namedListToTibble function returns the correct output NA in vector", {

    ##prepare input where only the list elements are named
    l <- list(a=LETTERS[1:10], b=c(letters[1:5], rep(NA, 5)))

    #setup expected data
    expected <- tibble(
      names = c(rep("a", 10), rep("b", 10)),
      variables = c(LETTERS[1:10], letters[1:5], rep(NA, 5))
    )

  ##run function
  output <- namedListToTibble(l)

  ##test
  expect_identical(expected, output)
})

test_that("check that the namedListToTibble function returns the correct output NA in list element names", {

    ##prepare input where only the list elements are named
    l <- list(LETTERS[1:10], letters[1:10])
    names(l) <- c("a")

    #setup expected data
    expected <- tibble(
      names = c(rep("a", 10), rep(NA, 10)),
      variables = c(LETTERS[1:10], letters[1:10])
    )

  ##run function
  output <- namedListToTibble(l)

  ##test
  expect_identical(expected, output)
})

test_that("check that the matrix_to_tibble function returns the correct output 1", {

  ####TEST1####
  ##prepare normal input data with rownames
  nd <- rnorm(20)
  m <- matrix(
        nd,
        ncol = 2,
        dimnames = list(letters[1:10],
        LETTERS[1:2])
    )

  #setup expected data
  expected <- tibble(
    my.rownames = letters[1:10],
    A = nd[1:10],
    B = nd[11:20]
  )

  ##run function
  output <- matrix_to_tibble(m, rowname = "my.rownames")

  ##test
  expect_equal(output, expected)
})

test_that("check that the matrix_to_tibble function returns the correct output 2", {

  ####TEST1####
  ##prepare normal input data without rownames
  nd <- rnorm(20)
  m <- matrix(nd, ncol = 2)
  colnames(m) <- LETTERS[1:2]
  
  #setup expected data
  expected <- tibble(
    A = nd[1:10],
    B = nd[11:20]
  )

  ##run function
  output <- matrix_to_tibble(m)

  ##test
  expect_equal(output, expected)
})

test_that("check that the matrix_to_tibble function returns the correct output 3", {

  ####TEST1####
  ##prepare normal input data without dim names
  nd <- rnorm(20)
  m <- matrix(nd, ncol = 2)
  
  #setup expected data
  expected <- tibble(
    V1 = nd[1:10],
    V2 = nd[11:20]
  )

  ##run function
  output <- matrix_to_tibble(m)

  ##test
  expect_equal(output, expected)
})
