alphaDiag <-
function(X, library = "GUDHI", printProgress = FALSE) {

  if (!is.numeric(X) && !is.data.frame(X)) {
    stop("X should be a matrix of coordinates")
  }
  if (library == "gudhi" || library == "Gudhi") {
    library <- "GUDHI"
  }
  if (library != "GUDHI") {
    stop("library should be 'GUDHI'")
  }
  if (!is.logical(printProgress)) {
    stop("printProgress should be logical")
  }

	if (library=="GUDHI")
	{
		max_num_pairs=5000  # to be added as an option
		points=as.vector(t(X))
		dim=ncol(X)
		num_points=nrow(X)
		diagram=matrix(-1, ncol=3, nrow=max_num_pairs)
		
		out=.C("alphashape_persistence_diagram_GUDHI", as.double(points),
	                                       as.integer(dim),
	                                       as.integer(num_points),
	                                       NCOL(X) - 1,
	                                       as.double(diagram),
	                                       as.integer(printProgress),
	                                       dup=FALSE, package="TDA")
		## GUDHI in windows has a problem: instead of Inf writes 1.#INF
		## here we manually fix the problem
		Diag=read.csv("outputTDA.txt", sep=" ", header=F)
		attributes(Diag[,3])$levels=c(attributes(Diag[,3])$levels, "Inf")
		change=which(Diag[,3]=="1.#INF")
		Diag[change,3]="Inf"
		Diag[,3]=as.numeric(as.character(Diag[,3]))
		Diag=as.matrix(Diag)
	}
	
		N=dim(Diag)[1]
		remove=NULL  # we remove points with lifetime=0
		for (i in 1:N){
			if (Diag[i,2]==Diag[i,3]) remove=c(remove,i)
		}
		#remove points with lifetime=0
		if (!is.null(remove)) Diag=Diag[-remove,]  
		## change Inf values to maxscale
		#Diag[which(Diag[,3]==Inf),3]=maxscale	

		colnames(Diag)=c("dimension","Birth", "Death")
		
		if (class(Diag)!="matrix") Diag=t(Diag) #in the case there is only 1 point
		
		class(Diag)="diagram"
		#attributes(Diag)$maxdimension=maxdimension
		#attributes(Diag)$scale=c(0, maxscale)
  attributes(Diag)[["maxdimension"]] <- max(Diag[, 1])
  nonInf <- which(Diag[, 2] != Inf & Diag[, 3] != Inf)
  attributes(Diag)[["scale"]] <-
      c(min(Diag[nonInf, 2:3]), max(Diag[nonInf, 2:3]))
		attributes(Diag)$call=match.call()
			
		return(list("diagram"=Diag))
	
}
