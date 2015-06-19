\name{alphaDiag}
\alias{alphaDiag}
\title{
Persistence Diagram of Alpha Shape
}
\description{
\code{alphaDiag} function computes the persistence diagram of the alpha shape filtration built on top of a point cloud.
}
\usage{
alphaDiag(X, maxdimension = NCOL(X) - 1, maxscale, library = "GUDHI",
          printProgress = FALSE)
}
\arguments{
  \item{X}{
    an \eqn{n} by \eqn{d} matrix of coordinates, used by the function \code{FUN}, where \eqn{n} is the number of points stored in \code{X} and \eqn{d} is the dimension of the space.
}
  \item{maxdimension}{
    a number that indicates the maximum dimension of the homological features to compute: 0 for connected components, 1 for loops, 2 for voids, and so on. Default is \eqn{d-1}, which is (dimension of embedding space - 1).
}
  \item{maxscale}{
    number: maximum value of the alpha shape filtration.
}
  \item{library}{
    Currently only the library 'GUDHI' can be used to compute the persistence diagram.
}
  \item{printProgress}{
	  if \code{TRUE}, a progress bar is printed. Default is \code{FALSE}.
}
}
\details{
\code{alphaDiag} constructs the alpha shape filtration, using the C++ library GUDHI.
See refereneces.
}
\value{
\code{alphaDiag} returns a list with the following elements:
  \item{diagram}{
    an object of class \code{diagram}, a \eqn{P} by 3 matrix, where \eqn{P} is the number of points in the resulting persistence diagram. The first column stores the dimension of each feature (0 for components, 1 for loops, 2 for voids, etc). Second and third columns are Birth and Death of the features.
}
}
\references{
Kaspar Fischer, Introduction to Alpha Shapes, http://www.cs.uu.nl/docs/vakken/ddm/texts/Delaunay/alphashapes.pdf

Jean-Daniel Boissonnat, Marc Glisse, Clement Maria, and Vincent Rouvreau, "GUDHI, Simplicial Complexes and Persistent Homology Packages". https://project.inria.fr/gudhi/software/.

Dmitriy Morozov, Homological Illusions of Persistence and Stability, http://www.mrzv.org/publications/thesis/phd/

Herbert Edelsbrunner and John Harer (2010), Computational topology: an introduction. American Mathematical Society.
}
\author{
Jisu Kim and Vincent Rouvreau
}
\seealso{
\code{\link{summary.diagram}}, \code{\link{plot.diagram}}, \code{\link{gridDiag}}, \code{\link{ripsDiag}}
}
\examples{
}
\keyword{methods}

