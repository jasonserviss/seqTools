context("generalTools methods")


##run test
test_that("check that the namedListToTibble function returns the correct output", {
  
  ####TEST1####
  ##prepare normal input data
  l <- list(a=LETTERS[1:10], b=letters[1:5])
  
  #setup expected data
  expected <- tibble(
    names
  )
  
  ##run function
  codedMat <- cpm(data)
  
  ##test
  expect_true()
  expect_equal()
  expect_equal()
})


##run test
test_that("check that the matrix_to_tibble function returns the correct output", {
  
  ####TEST1####
  ##prepare normal input data
  m <- matrix(
        rnorm(20),
        ncol = 2,
        dimnames = list(letters[1:10],
        LETTERS[1:2])
    )
  
  #setup expected data
  expected <- tibble(
    names
  )
  
  ##run function
  output <- matrix_to_tibble(m)
  
  ##test
  expect_true()
  expect_equal()
  expect_equal()
})
