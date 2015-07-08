ripsDiag <-
function(X, maxdimension, maxscale, dist = "euclidean", library = "GUDHI",
         printProgress = FALSE) {

  if (!is.numeric(X) && !is.data.frame(X)) {
    stop("X should be a matrix of coordinates")
  }
  tryCatch(maxdimension <- as.double(maxdimension), error = function(e) {
      stop("maxdimension should be numeric")})
  if (length(maxdimension) != 1 || maxdimension < 0) {
    stop("maxdimnsion should be a nonnegative integer")
  }
  tryCatch(maxscale <- as.double(maxscale), error = function(e) {
      stop("maxscale should be numeric")})
  if (length(maxscale) != 1) {
    stop("maxscale should be a number")
  }
  if (dist != "euclidean" && dist != "arbitrary") {
    stop ("dist should be either 'euclidean' or 'arbitrary'")
  }
  if (library == "dionysus" || library == "DIONYSUS") {
    library <- "Dionysus"
  }
  if (library == "gudhi" || library == "Gudhi") {
    library <- "GUDHI"
  }
  if (library != "Dionysus" && library != "GUDHI") {
    stop("library should be either 'Dionysus' or 'GUDHI'")
  }
  if (dist == "arbitrary") {
    library <- "Dionysus"
  }
  if (!is.logical(printProgress)) {
    stop("printProgress should be logical")
  }

  X <- as.matrix(X)

  if (dist == "arbitrary") {
    library <- "Dionysus"
  }
  # in 32bit architectures Dionysus L2 doesn't work
  if (dist == "euclidean" && library == "Dionysus") {
    dist <- "arbitrary"
    X <- as.matrix(dist(X))
  }

  max_num_pairs <- 5000  # to be added as an option

  ripsOut <- RipsDiag(X = X, maxdimension = maxdimension, maxscale = maxscale,
      dist = dist, library = library, printProgress = printProgress)

  Diag <- ripsOut[[1]]
  if (NROW(Diag) > 0) {
    ## change Inf values to maxscale
    Diag[which(Diag[, 3] == Inf), 3] <- maxscale  
  }

  colnames(Diag) <- c("dimension", "Birth", "Death")
  class(Diag) <- "diagram"
  attributes(Diag)[["maxdimension"]] <- max(Diag[, 1])
  attributes(Diag)[["scale"]] <- c(0, maxscale)
  attributes(Diag)$call <- match.call()
  out <- list("diagram" = Diag)
  return (out)
}
