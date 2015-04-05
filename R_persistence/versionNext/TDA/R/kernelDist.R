kernelDist <-
function(X, Grid, h, printProgress = FALSE) {

  if (!is.numeric(X) && !is.data.frame(X)) {
    stop("X should be a matrix of coordinates")
  }
  if (!is.numeric(Grid) && !is.data.frame(Grid)) {
    stop("Grid should be a matrix of coordinates")
  }
  if (!is.vector(h) || length(h) != 1) {
    stop("h should be a positive number")
  }
  if (!is.logical(printProgress)) {
    stop("printProgress should be a logical variable")
  }
  
  return (KdeDist(X = as.matrix(X), Grid = as.matrix(Grid), h = as.double(h), printProgress = printProgress))
}
