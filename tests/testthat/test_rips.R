context("ripsFiltration")

test_that("default circle example ripsFiltration", {
    n <- 5
    X <- cbind(cos(2*pi*seq_len(n)/n), sin(2*pi*seq_len(n)/n))
    maxdimension <- 1
    maxscale <- 1.5
    FltRips <- ripsFiltration(X = X, maxdimension = maxdimension,
               maxscale = maxscale, dist = "euclidean", library = "D2",
               printProgress = TRUE)
    expect_equal(FltRips$cmplx[[1]],1)
    expect_true(FltRips$increasing)
    expect_equal(FltRips$values[[1]],0)
    expect_true(abs(FltRips$values[[8]]-1.175571)<.000001)
})

test_that("One dimensional ripsFiltration in a line", {
    Y  <-  matrix(c(1,2.1,3.3,4.6,6))
    FltRips  <- ripsFiltration(X =Y, maxdimension = 1,
                 maxscale = 1.5, dist = "euclidean", library = "D2",
                 printProgress = TRUE)
    expect_equal(FltRips$cmplx[[1]],1)
    expect_equal(FltRips$cmplx[[6]],c(1,2))
    expect_equal(FltRips$cmplx[[9]],c(4,5))
    expect_equal(FltRips$cmplx[[8]],c(3,4))
    expect_true(FltRips$increasing)
})


