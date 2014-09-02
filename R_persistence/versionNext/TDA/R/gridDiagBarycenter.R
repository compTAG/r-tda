gridDiagBarycenter <-
function(X, FUN, lim, by=(lim[2]-lim[1])/20, sublevel=TRUE, printStatus=FALSE, diagLimit=NULL, ...){

	if (!is.function(FUN)) stop("FUN should be a function")	
	if (!is.numeric(X) && !is.data.frame(X)) stop("X should be a matrix of coordinates")
  if (2*ncol(X)!=length(lim)) stop("dimension of X does not match with lim")
	if (!is.vector(by) || length(by)!=1) stop("by should be a positive number")
	if (!is.logical(sublevel)) stop("sublevel should be logical")
	if (!is.logical(printStatus)) stop("printStatus should be logical")
	if (!is.null(diagLimit) && (!is.vector(diagLimit) || length(diagLimit)!=1) ) stop("diagLimit should be a positive number")	

	# in case there is only 1 point
	if (is.vector(X)) X=t(X)
	
#	if (ncol(X)>3) 

	Grid=gridByBarycenter(lim, by=by)
	p=FUN(X,Grid$grid,...)
		
	dim=Grid$dim

  gridValues=p
  if (sublevel==FALSE) gridValues=max(p)-p
	
	#write input.txt and read output.txt
  if (ncol(X)<=3)
  {
  	computeGrid=.C("gridMem",extFcnVal=as.double(gridValues),extDim=as.integer(length(dim)),extGridNum=as.integer(dim),input=as.integer(printStatus),
                   dup=TRUE, package="persistence")
  }
  else
  {
    computeGrid=.C("gridBarycenter",extFcnVal=as.double(gridValues),extDim=as.integer(length(dim)),extGridNum=as.integer(dim),input=as.integer(printStatus),
                   dup=TRUE, package="persistence")    
  }
	out=read.table("outputDionysus.txt", sep="\n")
	
	##convert output.txt in matrix format
	vecOut=as.vector(out$V1)
	whichDimens=c((1:length(vecOut))[-grep(" ",vecOut)], length(vecOut)+1)
	dim=NULL
	for (i in 1:(length(whichDimens)-1)){
		dim=c(dim, rep(i-1, whichDimens[i+1]-whichDimens[i]-1))
		}
	out2=data.frame(dim,life=out[grep(" ",vecOut),1])
	life2=matrix(NA, ncol=2, nrow=length(dim))
	for (i in 1:length(dim)){
		life2[i,]=as.numeric(unlist(strsplit(as.character(out2[i,2]), " "))	)
		}
	
	Diag=cbind(dim,life2)
	if (sublevel==FALSE) {
		colnames(Diag)=c("dim", "Death", "Birth")
		Diag[,2:3]=max(p)-Diag[,2:3]
		Diag[1,3]= min(p) #check
		Diag[1,2]= ifelse(is.null(diagLimit), Diag[1,2], diagLimit) #check
	} else {
		colnames(Diag)=c("dim", "Birth", "Death")
		Diag[1,3]=ifelse(is.null(diagLimit), max(p), diagLimit) 
	}
	if (sublevel==FALSE) Diag[,2:3]=Diag[,3:2]

	if (class(Diag)!="matrix") Diag=t(Diag) #in the case there is only 1 point
	class(Diag)="diagram"
	attributes(Diag)$maxdimension=max(Diag[,1])
	nonInf=which(Diag[,2]!=Inf & Diag[,3]!=Inf)
	attributes(Diag)$scale=c(min(Diag[nonInf,2:3]), max(Diag[nonInf,2:3]))
	attributes(Diag)$call=match.call()
	return(Diag)
}
