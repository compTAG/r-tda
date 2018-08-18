context("FiltrationDiag")

test_that("FiltrationDiag works with dionysus2" , {
    X <- matrix(c(0,0,100,100,0,102,0,101),nrow=4)
    Fltrips = ripsFiltration(X,maxdimension = 1, maxscale = 120, library = "Dionysus")
    DiagRips = filtrationDiag(Fltrips, maxdimension = 0, library = "Dionysus")
    DiagRips2 = filtrationDiag(Fltrips, maxdimension = 0, library = "D2", location = FALSE)
    for (i in 1:nrow(DiagRips)) {
        for (j in 1:ncol(DiagRips)) {
            expect_equal(DiagRips[i,j],DiagRips2[i,j])
        }
    }
})


