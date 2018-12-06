#' @param X an n by d matrix of coordinates containing points used int the
#' density estimation, where n is the number of points and d is the dimension
#' @param Grid an m by d matrix of coordinates, where m is the number of points
#' in the grid and d is the dimension
#' @param k a number used as the smoothing parameter
#' @return a vector of length m (number of points in Grid) containing the value
#' of the knn Density Estimator for each point in the grid

knnDE <-
function(X, Grid, k){
  
  # ensure input X is a matrix of numbers
  if (!is.numeric(X) && !is.data.frame(X)) {
    stop("X should be a matrix of coordinates")
  }
  # ensure input Grid is a matrix of numbers
  if (!is.numeric(Grid) && !is.data.frame(Grid)) {
    stop("Grid should be a matrix of coordinates")
  }
  # ensure that k is a positive integer
  if (!is.numeric(k) || length(k) != 1 || k < 1) {
    stop("k should be a positive integer")
  }

  # store the number of rows and columns in input matrix X
  # (d being the dimension and n being number of points in X)
  d <- ncol(X)
  n <- nrow(X)
  # find the Euclidean distances from k nearest neighbors using input matrix X, the query matrix
  # Grid, a maximum number n of nearest neighbors to search, and the nearest neighbor searching
  # algorithm ("kd_tree")
  r.k <- apply(FNN::knnx.dist(X, Grid, k = k, algorithm = "kd_tree"), 1, max)
  # volume of the Euclidean d dimensional unit ball
  v.d <- pi^(d/2) /gamma(d/2+1)
  # function in its final form, dividing smoothing parameter k by the number of points in X multiplied 
  # by the volume of unit ball and the Euclidean distances from k nearest neighbors to the power of d dimensions
  out <- k / (n * v.d * r.k ^ d)  
  return(out)
}
