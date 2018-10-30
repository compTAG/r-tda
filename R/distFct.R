#'Distance Function
#'
#'@description The function distFct computes the distance between each point of a set
#' Grid and the corresponding closest point of another set X.
#'
#'@param X a numeric m by d matrix of coordinates in the space, where m is the number
#'of points in X and d is the dimension of the space. X is the set of points whose distance
#'is being measured from a base Grid.
#'
#'@param Grid a numeric n by d matrix of coordinates in the space, where n is the number of
#'points in Grid and d is the dimension of the space. Grid serves as a base set where each point in
#'Grid is compared to the closest point in X.
#'
#' @details Given a set of points X, the distance function computed at g is defined as
#'
#'  \deqn{d(g) = \inf_{x\in X} ||x-g||_2}
#'
#'@return The function distFct returns a numeric vector V of length n, where n is the number
#'  of points stored in Grid. Each value in V corresponds to the distance between a point in G
#'  and the nearest point in X.
#'
#'@author Fabrizio Lecci
#'
#'@seealso kde,kernelDist,dtm
#'
#'@examples
#'  ## Generate Data from the unit circle
#'  n <- 300
#'  X <- circleUnif(n)
#'
#'  ## Construct a grid of points over which we evaluate the function
#'  by <- 0.065
#'  Xseq <- seq(-1.6, 1.6, by = by)
#'  Yseq <- seq(-1.7, 1.7, by = by)
#'  Grid <- expand.grid(Xseq, Yseq)
#'
#'  ## distancefct
#'  distance <- distFct(X, Grid)

distFct <-
function(X, Grid) {

  if (!is.numeric(X) && !is.data.frame(X)) {
    stop("X should be a matrix of coordinates")
  }
  if (!is.numeric(Grid) && !is.data.frame(Grid)) {
    stop("Grid should be a matrix of coordinates")
  }
  if (NCOL(X) != NCOL(Grid)) {
    stop("dimension of X does not match with dimension of Grid")
  }
    
  X <- as.matrix(X) 
  Grid <- as.matrix(Grid)
  distances <- FNN::knnx.dist(X, Grid, k = 1, algorithm = c("kd_tree")) 
  dOut <- as.vector(distances)
  return(dOut)
}
