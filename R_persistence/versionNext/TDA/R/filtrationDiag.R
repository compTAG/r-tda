filtrationDiag <- function(
    filtration, maxdimension = max(sapply(X = filtration[[1]], FUN = length, simplify = TRUE)) - 1, library = "Dionysus", location = FALSE, printProgress = FALSE, diagLimit = NULL) {

  filtrationOut <- FiltrationDiag(
      filtration = filtration, maxdimension = maxdimension, library = library,
	  location = location, printProgress = printProgress)

   if (location == TRUE) {
    BirthLocation <- filtrationOut[[2]][, 1]
    DeathLocation <- filtrationOut[[2]][, 2]
    if (library == "Dionysus") {
      CycleLocation <- filtrationOut[[3]]
    }
  }

  Diag <- filtrationOut[[1]]
  if (NROW(Diag) > 0) {
    Diag[Diag == Inf] <- ifelse(is.null(diagLimit), Inf, diagLimit) 
  }

  colnames(Diag) <- c("dimension", "Birth", "Death")

  class(Diag) <- "diagram"
  attributes(Diag)[["maxdimension"]] <- max(Diag[, 1])
  nonInf <- which(Diag[, 2] != Inf & Diag[, 3] != Inf)
  attributes(Diag)[["scale"]] <-
    c(min(Diag[nonInf, 2:3]), max(Diag[nonInf, 2:3]))  
  attributes(Diag)[["call"]] <- match.call()
  if (location == FALSE || library == "GUDHI") {
    out <- list("diagram" = Diag)
  } else if (library == "PHAT") {
    out <- list(
	    "diagram" = Diag, "birthLocation" = BirthLocation,
		"deathLocation" = DeathLocation)
  } else {
    out <- list(
	    "diagram" = Diag, "birthLocation" = BirthLocation,
		"deathLocation" = DeathLocation, "cycleLocation" = CycleLocation)
  }
  return (out)
}