landscape <-
function(Diag, dimension = 1, KK = 1,
         tseq = seq(min(Diag[, 2:3]), max(Diag[, 2:3]), length = 500)) {
    
  if (class(Diag) != "diagram" && class(Diag) != "matrix" &&
      !is.data.frame(Diag)) {
    stop("Diag should be a diagram, or a P by 3 matrix")
  }
  tryCatch(dimension <- as.double(dimension), error = function(e) {
      stop("dimension should be numeric")})
  if (length(dimension) != 1 || dimension < 0) {
    stop("dimension should be an nonnegative integer")
  }
  tryCatch(KK <- as.double(KK), error = function(e) {
      stop("KK should be numeric")})
  if (min(KK) <= 0) {
    stop("KK should be positive integer")
  }
  tryCatch(tseq <- as.double(tseq), error = function(e) {
      stop("tseq should be numeric")})
    
    isNA=length(which(Diag[, 1] == dimension))
    if (isNA==0) return(rep(0, length(tseq))) #in case there are no features with dimension "dimension"
    	
	Diag=Diag[which(Diag[,1]==dimension),]
	if (class(Diag)!="matrix") Diag=t(Diag) #in the case there is only 1 point
	
	Npoints=nrow(Diag)

    fab = matrix(NA, nrow = length(tseq), ncol = Npoints)
    lambda = numeric()
    for (j in 1:Npoints) {    
        fab[,j]=sapply(1:length(tseq), FUN=function(i){
        	max(min(tseq[i] - Diag[j, 2], Diag[j,3] - tseq[i]), 0)        	    
        })        
    }
    lambda=sapply(1:length(tseq),  FUN=function(i){
    	sort(fab[i, ], decreasing = TRUE)[KK]  	
    })
  lambda[is.na(lambda)] <- 0
  if (length(KK) == 1) {
    lambda <- matrix(lambda)
  } else {
    lambda <- t(lambda)
  }
    return(lambda)
}