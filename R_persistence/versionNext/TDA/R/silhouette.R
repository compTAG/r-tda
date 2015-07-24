silhouette <-
function(Diag, p = 1, dimension = 1,
         tseq = seq(min(Diag[, 2:3]), max(Diag[, 2:3]), length = 500)) {

  if ((class(Diag) != "diagram" && class(Diag) != "matrix" &&
      !is.data.frame(Diag)) || NCOL(Diag) != 3) {
    stop("Diag should be a diagram, or a P by 3 matrix")
  }
  if (!is.numeric(dimension) || length(dimension) != 1 || dimension < 0) {
    stop("dimension should be an nonnegative integer")
  }
  if (!is.numeric(p) || length(p) != 1 || any(p < 0)) {
    stop("p should be a nonnegative number")
  }
	if (!is.numeric(tseq)) {
    stop("tseq should be a numeric vector")
  }
  
  isNA <- length(which(Diag[, 1] == dimension))
  if (isNA == 0) {
    return(rep(0, length(tseq))) #in case there are no features with dimension "dimension"
  }
    	
	Diag <- Diag[which(Diag[,1] == dimension), ]
	if (class(Diag) != "matrix") {
    Diag <- t(Diag) #in the case there is only 1 point
  }

	left <- Diag[, 2]
	right <- Diag[, 3]
	Npoints <- length(left)
   
	### Silhouette
	w <- (right - left) ^ p
	w <- w / sum(w)
   
	ff <- rep(0, length(tseq))

	for(i in seq_len(Npoints)) {
		tmp <- pmax(pmin(tseq - left[i], right[i] - tseq), 0)
		ff <- ff + tmp * w[i]
	}
	return(ff)
}
