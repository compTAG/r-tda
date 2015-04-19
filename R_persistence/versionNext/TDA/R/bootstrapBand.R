bootstrapBand <- 
function(X, FUN, Grid, B = 30, alpha = 0.05, parallel = FALSE,
         printProgress = FALSE, weight = NULL, ...) {
     
  if (!is.numeric(X) && !is.data.frame(X)) {
    stop("X should be a matrix of coordinates")
  }
  if (class(FUN) != "function") {
    stop("FUN should be function")
  }
  if (!is.numeric(Grid) && !is.data.frame(Grid)) {
    stop("Grid should be a matrix of coordinates")
  }
  tryCatch(B <- as.double(B), error = function(e) {
      stop("B should be numeric")})
  if (length(B) != 1 || B < 1) {
    stop("B should be a positive integer")
  }
  tryCatch(alpha <- as.double(alpha), error = function(e) {
      stop("alpha should be numeric")})
  if (alpha < 0 || alpha > 1) {
    stop("alpha should be a number between 0 and 1")
  }
  if (!is.logical(parallel)) {
    stop("parallel should be logical")
  }
  if (!is.logical(printProgress)) {
    stop("printProgress should be logical")
  }
  if (((length(weight) != 1 && length(weight) != NROW(X)) ||
      !is.numeric(weight)) && !is.null(weight)) {
    stop("weight should be either NULL, a number, or a vector of length equals the number of sample")
  }

  X <- as.matrix(X)

  if (is.null(weight)) {
    ff <- FUN(X, Grid, ...)
    boostFUN <- function(i) {
      I <- sample(NROW(X), replace = TRUE, size = NROW(X))
      bootF <- FUN(X[I, , drop = FALSE], Grid, ...)
      width1 <- max(abs(ff - bootF))
      if (printProgress) {
        cat(i," ")
      }
      return (width1)
    }
  } else {
    ff <- FUN(X, Grid, weight = weight, ...)
    boostFUN <- function(i) {
      weightBoost <- rMultinom(size = sum(weight), prob = weight)
      bootF <- FUN(X[I, , drop = FALSE], Grid, weight = weightBoost, ...)
      width1 <- max(abs(ff - bootF))
      if (printProgress) {
        cat(i," ")
      }
      return (width1)
    }
  }
  if (parallel) {
    boostLapply <- mclapply
  } else {
    boostLapply <- lapply
  }

  if (printProgress) {
    cat("Bootstrap: ")
  }
  width <- boostLapply(seq_len(B), FUN = boostFUN)
  if (printProgress) {
    cat("\n")
  }
  width <- quantile(unlist(width), 1 - alpha)
     
  UPband <- ff + width
  LOWband <- ff - width
  # LOWband[which(LOWband < 0)] = 0   #set negative values of lower band =0
  Band <- cbind(LOWband, UPband)

  return (list("width" = width, "fun" = ff, "band" = Band))
}