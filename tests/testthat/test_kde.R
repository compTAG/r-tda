context("KDE")

test_that("kde returns 1 when points are the same", {
    h <- 1/sqrt(2*pi)
    X  <-  expand.grid(0,0);Grid  <-  expand.grid(0,0)
    expect_equal(kde(X=X,Grid=Grid,h=1/sqrt(2*pi)),1)
})

test_that("peaks are higher than surrounding points", {
    h <- 1/sqrt(2*pi)
    Grid <- expand.grid(seq(-2,2,by = .01)) 
    X <- expand.grid(c(-1,1))
    KDE = kde(X = X, Grid = Grid, h=h)
    checks = c(max(KDE) > KDE[1], max(KDE) > KDE[401], max(KDE) > KDE[201])
    for (bool in checks) {
        expect_true(bool)
    }
})


