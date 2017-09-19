funFiltration <- function(FUNvalues, cmplx) {

  funOut <- FunFiltration(FUNvalues = FUNvalues, cmplx = cmplx)
  
  out <- list("cmplx" = funOut[[1]], "values" = funOut[[2]])
  return (out)
}