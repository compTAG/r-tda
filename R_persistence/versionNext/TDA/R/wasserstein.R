wasserstein <-
function(Diag1, Diag2, p = 1, dimension = 1){

  if (class(Diag1) != "diagram" && class(Diag1) != "matrix" && !is.data.frame(Diag1))
    stop("Diag1 should be a diagram or a matrix")
  if (class(Diag2) != "diagram" && class(Diag2) != "matrix" && !is.data.frame(Diag1))
    stop("Diag2 should be a diagram or a matrix")
  if (!is.vector(p) || length(p) != 1 || p < 1)
    stop("p should be a positive integer")
  if (!is.vector(dimension) || !all(dimension >= 0))
    stop("dimension should be a nonnegative integer or a vector of nonnegative integer")

  wassersteinDistance <- rep(0, length(dimension))
  for (dimIdx in seq(along = dimension)) {
    Diag1Dim <- Diag1[which(Diag1[, 1] == dimension[dimIdx]), 2:3, drop = FALSE]
    Diag2Dim <- Diag2[which(Diag2[, 1] == dimension[dimIdx]), 2:3, drop = FALSE]

    if (length(Diag1Dim) == 0) {
      Diag1Dim <- array(c(dimension[dimIdx], 0, 0), dim = c(1, 3))
    }
    if (length(Diag2Dim) == 0) {
      Diag2Dim <- array(c(dimension[dimIdx], 0, 0), dim = c(1, 3))
    }

    wassersteinDistance[dimIdx] <- Wasserstein(Diag1Dim, Diag2Dim, p)
  }

  return (max(wassersteinDistance))
}
