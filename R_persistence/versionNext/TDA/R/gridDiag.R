gridDiag <-
function(X, FUN, lim, by, maxdimension = length(lim) / 2 - 1, sublevel = TRUE,
         library = "Dionysus", location = FALSE, printProgress = FALSE,
         diagLimit = NULL, ...) {

  if (!is.numeric(X) && !is.data.frame(X)) {
    stop("X should be a matrix of coordinates")
  }
  if (!is.function(FUN)) {
    stop("FUN should be a function")
  }
  tryCatch(lim <- as.double(lim), error = function(e) {
      stop("lim should be numeric")})
  if (length(lim) %% 2 != 0) {
    stop("lim should be either a matrix or a vector of even elements")
  }
  if (2 * NCOL(X) != length(lim)) {
    stop("dimension of X does not match with lim")
  }
  tryCatch(by <- as.double(by), error = function(e) {
      stop("by should be numeric")})
  if ((length(by) != 1 && length(by) != NCOL(X)) || min(by) <= 0) {
    stop("by should be either a positive number or a positive vector of length equals dimension of grid")
  }
  tryCatch(maxdimension <- as.double(maxdimension), error = function(e) {
      stop("maxdimension should be numeric")})
  if (length(maxdimension) != 1 || maxdimension < 0) {
    stop("maxdimnsion should be a nonnegative integer")
  }
  if (!is.logical(sublevel)) {
    stop("sublevel should be logical")
  }
  if (library == "dionysus" || library == "DIONYSUS") {
    library <- "Dionysus"
  }
  if (library == "phat" || library == "Phat") {
    library <- "PHAT"
  }
  if (library != "Dionysus" && library != "PHAT") {
    stop("library should be a string: either 'Dionysus' or 'PHAT'")
  }
  if (!is.logical(location)) {
    stop("location should be logical")
  }
  if (!is.logical(printProgress)) {
    stop("printProgress should be logical")
  }
  if ((!is.vector(diagLimit) || length(diagLimit) != 1) &&
      !is.null(diagLimit)) {
    stop("diagLimit should be a positive number")
  }

  X <- as.matrix(X)
  maxdimension <- min(maxdimension, NCOL(X) - 1)

  Grid <- gridBy(lim = lim, by = by)
  FUNvalues <- FUN(X, Grid[["grid"]], ...)
  gridDim <- Grid[["dim"]]    
  if (sublevel == FALSE) {
    FUNvalues <- -FUNvalues
  }

  # compute persistence diagram of function values over a grid
  if (ncol(X) <= 3) {
    gridOut <- GridDiag(FUNvalues = FUNvalues, gridDim = as.integer(gridDim),
        maxdimension = as.integer(maxdimension), decomposition = "5tetrahedra",
        library = library, location = location, printProgress = printProgress)
  } else {
    gridOut <- GridDiag(FUNvalues = FUNvalues, gridDim = as.integer(gridDim),
        maxdimension = as.integer(maxdimension), decomposition = "barycenter",
        library = library, location = location, printProgress = printProgress)
  }

  if (location == TRUE) {
    BirthLocation <- Grid[["grid"]][gridOut[[2]][, 1], ]
    DeathLocation <- Grid[["grid"]][gridOut[[2]][, 2], ]
    if (library == "Dionysus")
    {
      CycleLocation <- lapply(gridOut[[3]], function(c) {Grid[["grid"]][c, ]})
    }
  }

  Diag <- gridOut[[1]]
  if (NROW(Diag) > 0) {
    Diag[1, 3] <- ifelse(is.null(diagLimit), max(FUNvalues), diagLimit) 
  }
  if (sublevel == FALSE) {
    colnames(Diag) <- c("dim", "Death", "Birth")
    Diag[, 2:3] <- -Diag[, 3:2]
  } else {
    colnames(Diag) <- c("dim", "Birth", "Death")
  }

  class(Diag) <- "diagram"
  attributes(Diag)[["maxdimension"]] <- max(Diag[, 1])
  nonInf <- which(Diag[, 2] != Inf & Diag[, 3] != Inf)
  attributes(Diag)[["scale"]] <-
    c(min(Diag[nonInf, 2:3]), max(Diag[nonInf, 2:3]))
  attributes(Diag)$call <- match.call()
  if (location == FALSE)
  {
    out <- list("diagram" = Diag)
  } else if (library == "PHAT")
  {
    out <- list("diagram" = Diag, "birthLocation" = BirthLocation,
        "deathLocation" = DeathLocation)
  } else
  {
    out <- list("diagram" = Diag, "birthLocation" = BirthLocation,
        "deathLocation" = DeathLocation, "cycleLocation" = CycleLocation)
  }
  return (out)
	
#   out=read.table("outputTDA.txt", sep="\n")
# 
#   ##convert outputTDA.txt in matrix format
#  	vecOut=as.vector(out$V1)
#   if (location==FALSE)
#   {
#  	  whichDimens=c((1:length(vecOut))[!grepl(" ",vecOut)], length(vecOut)+1)
#   } else
#   {
#     whichDimens=(1:length(vecOut))[!grepl("[ C]",vecOut)]
#   }
#  	dim=NULL
#    if (length(whichDimens)>1)
#    {
#    	for (i in 1:(length(whichDimens)-1)){
#    		dim=c(dim, rep(as.numeric(vecOut[whichDimens[i]]), whichDimens[i+1]-whichDimens[i]-1))
#    		}
#    }
#   if (location==FALSE)
#   {
#  	  life=out[grep(" ",vecOut),1]
#   } else
#   {
#     life=(out[grep(" ",vecOut),1])[1:length(dim)]
#   }
#   life2=matrix(as.numeric(unlist(strsplit(as.character(life)," "))),ncol=2, nrow=length(dim), byrow=TRUE)
# 	Diag=cbind(dim,life2)
#   if (location==TRUE)
#   {
#     loc=out[(grep("L",vecOut)+1):(grep("C",vecOut)-1),1]
#     loc2=matrix(as.numeric(unlist(strsplit(as.character(loc)," "))),ncol=2, nrow=length(loc), byrow=TRUE)
#     BirthLocation=Grid$grid[loc2[,1],]
#     DeathLocation=Grid$grid[loc2[,2],]
#     
#     if (library=="Dionysus")
#     {
#       cycle=c("",as.character(out[(grep("C",vecOut)+1):nrow(out),1]))
#       CycleLocation = lapply(strsplit(as.character(cycle)," "),function(c){Grid$grid[as.numeric(c),]})
#     }
#   }
# 	  
#   if (nrow(Diag)>0) {
#   	if (sublevel==FALSE) {
#   		colnames(Diag)=c("dim", "Death", "Birth")
#   		Diag[,2:3]=-Diag[,2:3]
#   		Diag[1,3]= ifelse(is.null(diagLimit), min(-FUNvalues), diagLimit)
#   	} else {
#   		colnames(Diag)=c("dim", "Birth", "Death")
#   		Diag[1,3]=ifelse(is.null(diagLimit), max(FUNvalues), diagLimit) 
#   	}
#   	if (sublevel==FALSE) Diag[,2:3]=Diag[,3:2]
#   }
# 
# 	if (class(Diag)!="matrix") Diag=t(Diag) #in the case there is only 1 point
# 	class(Diag)="diagram"
# 	attributes(Diag)$maxdimensionension=max(Diag[,1])
# 	nonInf=which(Diag[,2]!=Inf & Diag[,3]!=Inf)
# 	attributes(Diag)$scale=c(min(Diag[nonInf,2:3]), max(Diag[nonInf,2:3]))
# 	attributes(Diag)$call=match.call()
#   if (location==FALSE)
#   {
#     out=list("diagram"=Diag)
#   } else if (library=="PHAT")
#   {
#     out=list("diagram"=Diag,"birthLocation"=BirthLocation,"deathLocation"=DeathLocation)
#   } else
#   {
#     out=list("diagram"=Diag,"birthLocation"=BirthLocation,"deathLocation"=DeathLocation,"cycleLocation"=CycleLocation)
#   }
}
