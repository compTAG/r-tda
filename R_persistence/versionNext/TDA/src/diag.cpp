// for R
#include <R.h>
#include <R_ext/Print.h>

// for Rcpp
#include <Rcpp.h>

// for Rips
#include <tdautils/ripsL2.h>
#include <tdautils/ripsArbit.h>

// for grid
#include <tdautils/gridUtils.h>

// for kernel density
#include <tdautils/kernelUtils.h>

// for changing formats and typecasting
#include <tdautils/typecastUtils.h>

//for GUDHI
#include <tdautils/gudhiUtils.h>

// for Dionysus
#include <tdautils/dionysusUtils.h>

// for phat
#include <tdautils/phatUtils.h>



// GridDiag by Brittany T. Fasy
// modified by Jisu Kim for
// arbitrary dimension & using memory as an input & setting maximum dimension.
/** \brief Interface for R code, construct the persistence diagram
  * of sublevel/superlevel sets of a function evaluated over a grid of points
  *
  * @param[out] Rcpp::List     A list
  * @param[in]  FUNvalues      A vector of length m1*...*md of function values over grid
  * @param[in]  gridDim        A vector (m1, ..., md),
  *                            where mi is number of grid points in ith dimension
  * @param[in]  maxdimension   Max dimension of the homological features to be computed.
  *                            This equals (maximal dimension of the Rips complex) - 1
  * @param[in]  decomposition  Either "5tetrahedra" or "barycenter"
  * @param[in]  library        Either "Dionysus" or "PHAT"
  * @param[in]  location       Are location of birth point, death point,
  *                            and representative cycles returned?
  * @param[in]  printProgress  Is progress printed?
  */
 // [[Rcpp::export]]
Rcpp::List
GridDiag(const Rcpp::NumericVector & FUNvalues
       , const Rcpp::IntegerVector & gridDim
       , const int                   maxdimension
       , const std::string         & decomposition
       , const std::string         & library
       , const bool                  location
       , const bool                  printProgress
	) {
#ifdef LOGGING
	//rlog::RLogInit(argc, argv);

	stdoutLog.subscribeTo(RLOG_CHANNEL("topology/persistence"));
	//stdoutLog.subscribeTo(RLOG_CHANNEL("topology/chain"));
	//stdoutLog.subscribeTo(RLOG_CHANNEL("topology/vineyard"));
#endif

	std::vector< std::vector< std::vector< double > > > persDgm;
	std::vector< std::vector< std::vector< unsigned > > > persLoc;
	std::vector< std::vector< std::vector< std::vector< unsigned > > > > persCycle;

  Fltr filtration;

	// Generate simplicial complex from function values and grid
	if (decomposition[0] == '5') {
    simplicesFromGrid(filtration, FUNvalues, gridDim, maxdimension + 1);
	}
	if (decomposition[0] == 'b') {
    simplicesFromGridBarycenter(
        filtration, FUNvalues, gridDim, maxdimension + 1);
	}
	if (printProgress) {
    Rprintf("# Generated complex of size: %d \n", filtration.size());
	}

	// Sort the simplices with respect to function values
	filtration.sort(Smplx::DataComparison());

	// Compute the persistence diagram of the complex
	if (library[0] == 'D') {
    FiltrationDiagDionysus(
        filtration, maxdimension, location, printProgress, persDgm, persLoc,
        persCycle);
	}
  else if (library[0] == 'G') {

  }
  else {
    std::vector< phat::column > cmplx;
    std::vector< double > values;
    phat::boundary_matrix< phat::vector_vector > boundary_matrix;
    filtrationDionysusToPhat< phat::column, phat::dimension >(
      filtration, cmplx, values, boundary_matrix);
    FiltrationDiagPhat(
      cmplx, values, boundary_matrix, maxdimension, location, printProgress,
      persDgm, persLoc, persCycle);
  }

	// Output persistent diagram
	return Rcpp::List::create(
		concatStlToRcpp< Rcpp::NumericMatrix >(persDgm, true, 3),
		concatStlToRcpp< Rcpp::NumericMatrix >(persLoc, false, 2),
		StlToRcppMatrixList< Rcpp::List, Rcpp::NumericMatrix >(persCycle));
}



// [[Rcpp::export]]
double
Bottleneck(const Rcpp::NumericMatrix & Diag1
         , const Rcpp::NumericMatrix & Diag2
	) {
	return bottleneck_distance(RcppToDionysus< PersistenceDiagram<> >(Diag1),
			RcppToDionysus< PersistenceDiagram<> >(Diag2));
}



// [[Rcpp::export]]
double
Wasserstein(const Rcpp::NumericMatrix & Diag1
          , const Rcpp::NumericMatrix & Diag2
          , const int                   p
	) {
	return wasserstein_distance(RcppToDionysus< PersistenceDiagram<> >(Diag1),
			RcppToDionysus< PersistenceDiagram<> >(Diag2), p);
}



// KDE function on a Grid
// [[Rcpp::export]]
Rcpp::NumericVector
Kde(const Rcpp::NumericMatrix & X
  , const Rcpp::NumericMatrix & Grid
  , const double                h
  , const Rcpp::NumericVector & weight
  , const bool                  printProgress
	) {
	const double pi = 3.141592653589793;
	const unsigned dimension = Grid.ncol();
	const unsigned gridNum = Grid.nrow();
	const double den = pow(h, (int)dimension) * pow(2 * pi, dimension / 2.0);
	Rcpp::NumericVector kdeValue;
	int counter = 0, percentageFloor = 0;
	int totalCount = gridNum;

	if (printProgress) {
		printProgressFrame(Rprintf);
	}

	if (dimension <= 1) {
		kdeValue = computeKernel< Rcpp::NumericVector >(
				X, Grid, h, weight, printProgress, Rprintf, counter, totalCount,
				percentageFloor);
	}
	else {
		kdeValue = computeGaussOuter< Rcpp::NumericVector >(
				X, Grid, h, weight, printProgress, Rprintf, counter, totalCount,
				percentageFloor);

	}

	for (unsigned gridIdx = 0; gridIdx < gridNum; ++gridIdx) {
		kdeValue[gridIdx] /= den;
	}

	if (printProgress) {
		Rprintf("\n");
	}

	return kdeValue;
}



// kernel Dist function on a Grid
// [[Rcpp::export]]
Rcpp::NumericVector
KdeDist(const Rcpp::NumericMatrix & X
      , const Rcpp::NumericMatrix & Grid
      , const double                h
      , const Rcpp::NumericVector & weight
      , const bool printProgress
	) {
	const unsigned sampleNum = X.nrow();
	const unsigned dimension = Grid.ncol();
	const unsigned gridNum = Grid.nrow();
	// first = sum K_h(X_i, X_j), second = K_h(x, x), third = sum K_h(x, X_i)
	std::vector< double > firstValue;
	const double second = 1.0;
	std::vector< double > thirdValue;
	double firstmean;
	Rcpp::NumericVector kdeDistValue(gridNum);
	int counter = 0, percentageFloor = 0;
	int totalCount = sampleNum + gridNum;

	if (printProgress) {
		printProgressFrame(Rprintf);
	}

	firstValue = computeKernel< std::vector< double > >(
			X, X, h, weight, printProgress, Rprintf, counter, totalCount,
			percentageFloor);

	if (dimension <= 1) {
		thirdValue = computeKernel< std::vector< double > >(
				X, Grid, h, weight, printProgress, Rprintf, counter, totalCount,
				percentageFloor);
	}
	else {
		thirdValue = computeGaussOuter< std::vector< double > >(
				X, Grid, h, weight, printProgress, Rprintf, counter, totalCount,
				percentageFloor);
	}

	if (weight.size() == 1) {
		firstmean = std::accumulate(firstValue.begin(), firstValue.end(), 0.0) / sampleNum;
	}
	else {
		firstmean = std::inner_product(
				firstValue.begin(), firstValue.end(), weight.begin(), 0.0) / 
				std::accumulate(weight.begin(), weight.end(), 0.0);
	}

	for (unsigned gridIdx = 0; gridIdx < gridNum; ++gridIdx) {
		kdeDistValue[gridIdx] = std::sqrt(firstmean + second - 2 * thirdValue[gridIdx]);
	}

	if (printProgress) {
		Rprintf("\n");
	}

	return kdeDistValue;
}



// distance to measure function on a Grid
// [[Rcpp::export]]
Rcpp::NumericVector
Dtm(const Rcpp::NumericMatrix & knnDistance
  , const double                weightBound
  , const double                r
	) {
	const unsigned gridNum = knnDistance.nrow();
  unsigned gridIdx, kIdx;
  double distanceTemp = 0.0;
	Rcpp::NumericVector dtmValue(gridNum, 0.0);
  unsigned weightSumTemp;

  if (r == 2.0) {
    for (gridIdx = 0; gridIdx < gridNum; ++gridIdx) {
      for (kIdx = 0, weightSumTemp = 0; (double)weightSumTemp < weightBound;
          ++kIdx) {
        distanceTemp = knnDistance[gridIdx + kIdx * gridNum];
        dtmValue[gridIdx] += distanceTemp * distanceTemp;
        ++weightSumTemp;
      }
      dtmValue[gridIdx] += distanceTemp * distanceTemp *
          (weightBound - (double)weightSumTemp);
      dtmValue[gridIdx] = std::sqrt(dtmValue[gridIdx] / weightBound);
    }
  }
  else if (r == 1.0) {
    for (gridIdx = 0; gridIdx < gridNum; ++gridIdx) {
      for (kIdx = 0, weightSumTemp = 0; (double)weightSumTemp < weightBound;
          ++kIdx) {
        distanceTemp = knnDistance[gridIdx + kIdx * gridNum];
        dtmValue[gridIdx] += distanceTemp;
        ++weightSumTemp;
      }
      dtmValue[gridIdx] += distanceTemp *
          (weightBound - (double)weightSumTemp);
      dtmValue[gridIdx] /= weightBound;
    }
  }
  else {
    for (gridIdx = 0; gridIdx < gridNum; ++gridIdx) {
      for (kIdx = 0, weightSumTemp = 0; (double)weightSumTemp < weightBound;
          ++kIdx) {
        distanceTemp = knnDistance[gridIdx + kIdx * gridNum];
        dtmValue[gridIdx] += std::pow(distanceTemp, r);
        ++weightSumTemp;
      }
      dtmValue[gridIdx] += std::pow(distanceTemp, r) *
          (weightBound - (double)weightSumTemp);
      dtmValue[gridIdx] = std::pow(dtmValue[gridIdx] / weightBound, 1 / r);
    }
  }

  return (dtmValue);
}



// distance to measure function on a Grid, with weight
// [[Rcpp::export]]
Rcpp::NumericVector
DtmWeight(const Rcpp::NumericMatrix & knnDistance
        , const double                weightBound
        , const double                r
        , const Rcpp::NumericMatrix & knnIndex
        , const Rcpp::NumericVector & weight
  ) {
  const unsigned gridNum = knnDistance.nrow();
  unsigned gridIdx, kIdx;
  double distanceTemp = 0.0;
  Rcpp::NumericVector dtmValue(gridNum, 0.0);
  double weightTemp, weightSumTemp;

  if (r == 2.0) {
    for (gridIdx = 0; gridIdx < gridNum; ++gridIdx) {
      for (kIdx = 0, weightSumTemp = 0.0; weightSumTemp < weightBound;
          ++kIdx) {
        distanceTemp = knnDistance[gridIdx + kIdx * gridNum];
        weightTemp = weight[knnIndex[gridIdx + kIdx * gridNum] - 1];
        dtmValue[gridIdx] += distanceTemp * distanceTemp * weightTemp;
        weightSumTemp += weightTemp;
      }
      dtmValue[gridIdx] += distanceTemp * distanceTemp *
          (weightBound - weightSumTemp);
      dtmValue[gridIdx] = std::sqrt(dtmValue[gridIdx] / weightBound);
    }
  }
  else if (r == 1.0) {
    for (gridIdx = 0; gridIdx < gridNum; ++gridIdx) {
      for (kIdx = 0, weightSumTemp = 0.0; weightSumTemp < weightBound;
          ++kIdx) {
        distanceTemp = knnDistance[gridIdx + kIdx * gridNum];
        weightTemp = weight[knnIndex[gridIdx + kIdx * gridNum] - 1];
        dtmValue[gridIdx] += distanceTemp * weightTemp;
        weightSumTemp += weightTemp;
      }
      dtmValue[gridIdx] += distanceTemp * (weightBound - weightSumTemp);
      dtmValue[gridIdx] /= weightBound;
    }
  }
  else {
    for (gridIdx = 0; gridIdx < gridNum; ++gridIdx) {
      for (kIdx = 0, weightSumTemp = 0.0; weightSumTemp < weightBound;
          ++kIdx) {
        distanceTemp = knnDistance[gridIdx + kIdx * gridNum];
        weightTemp = weight[knnIndex[gridIdx + kIdx * gridNum] - 1];
        dtmValue[gridIdx] += std::pow(distanceTemp, r) * weightTemp;
        weightSumTemp += weightTemp;
      }
      dtmValue[gridIdx] += std::pow(distanceTemp, r) *
          (weightBound - weightSumTemp);
      dtmValue[gridIdx] = std::pow(dtmValue[gridIdx] / weightBound, 1 / r);
    }
  }

  return (dtmValue);
}



template< typename RcppList, typename RealVector >
void filtrationSort(RcppList & cmplx, RealVector & values) {

  std::vector< std::pair< double, std::pair< unsigned, RealVector > > >
      vipairs(cmplx.size());
  typename RcppList::iterator iCmplx;
  typename RealVector::iterator iValue;
  unsigned idx;
  typename std::vector< std::pair< double, std::pair< unsigned, RealVector > > >
    ::iterator iPair;

  iCmplx = cmplx.begin();
  iValue = values.begin();
  idx = 0;
  for (iPair = vipairs.begin(); iPair != vipairs.end();
       ++iPair, ++iValue, ++iCmplx, ++idx) {
    *iPair = std::make_pair(*iValue, std::make_pair(idx, *iCmplx));
  }
  std::sort(vipairs.begin(), vipairs.end());

  iCmplx = cmplx.begin();
  iValue = values.begin();
  for (iPair = vipairs.begin(); iPair != vipairs.end();
       ++iPair, ++iCmplx, ++iValue) {
    *iCmplx = (iPair->second).second;
    *iValue = iPair->first;
  }
}



// FiltrationDiag
/** \brief Interface for R code, construct the persistence diagram from the
 *         filtration.
 *
 * @param[out] Rcpp::List     A list
 * @param[in]  filtration     The input filtration
 * @param[in]  maxdimension   Max dimension of the homological features to be
 *                            computed.
 * @param[in]  library        Either "GUDHI", "Dionysus", or "PHAT"
 * @param[in]  location       Are location of birth point, death point, and
 *                            representative cycles returned?
 * @param[in]  printProgress  Is progress printed?
 */
// [[Rcpp::export]]
Rcpp::List FiltrationDiag(
    const Rcpp::List  & filtration,
    const int           maxdimension,
    const std::string & library,
    const bool          location,
    const bool          printProgress
) {

  std::vector< std::vector< std::vector< double > > > persDgm;
  std::vector< std::vector< std::vector< unsigned > > > persLoc;
  std::vector< std::vector< std::vector< std::vector< unsigned > > > >
      persCycle;

  Rcpp::List filtrationTemp(filtration);
  Rcpp::NumericVector values(filtration[1]);

  if (!std::is_sorted(values.begin(), values.end())) {
    Rcpp::List cmplx(filtration[0]);
    Rcpp::List cmplxTemp(cmplx.begin(), cmplx.end());
    Rcpp::NumericVector valuesTemp(values.begin(), values.end());
    filtrationSort(cmplxTemp, valuesTemp);
    filtrationTemp = Rcpp::List::create(cmplxTemp, valuesTemp);
  }

  if (library[0] == 'G') {
    int coeff_field_characteristic = 2;
    double min_persistence = 0.0;
    Gudhi::Simplex_tree<> smplxTree = filtrationRcppToGudhi<
        Gudhi::Simplex_tree<>, Rcpp::NumericVector >(filtrationTemp);
    FiltrationDiagGudhi(
        smplxTree, coeff_field_characteristic, min_persistence, maxdimension,
        printProgress, persDgm);
  }
  else if (library[0] == 'D') {
    FiltrationDiagDionysus(
      filtrationRcppToDionysus< Fltr, Rcpp::NumericVector >(filtrationTemp),
      maxdimension, location, printProgress, persDgm, persLoc, persCycle);
  }
  else {
    std::vector< phat::column > cmplx;
    std::vector< double > values;
    phat::boundary_matrix< phat::vector_vector > boundary_matrix;
    filtrationDionysusToPhat< phat::column, phat::dimension >(
      filtrationRcppToDionysus< Fltr, Rcpp::NumericVector >(filtrationTemp),
      cmplx, values, boundary_matrix);
    FiltrationDiagPhat(cmplx, values, boundary_matrix, maxdimension, location,
      printProgress, persDgm, persLoc, persCycle);
  }

  // Output persistent diagram
  return Rcpp::List::create(
    concatStlToRcpp< Rcpp::NumericMatrix >(persDgm, true, 3),
    concatStlToRcpp< Rcpp::NumericMatrix >(persLoc, false, 2),
    StlToRcppMatrixList< Rcpp::List, Rcpp::NumericMatrix >(persCycle));
}



// funFiltration
/** \brief Interface for R code, construct the persistence diagram from the
 *         filtration of the function values.
 *
 * @param[out] Rcpp::List  A list
 * @param[in]  FUNvalues   The inlut function values
 * @param[in]  cmplx       The input simplicial complex
 */
// [[Rcpp::export]]
Rcpp::List FunFiltration(
    const Rcpp::NumericVector & FUNvalues,
    const Rcpp::List          & cmplx
) {

  const unsigned nCmplx = cmplx.size();
  Rcpp::List rcppCmplx(cmplx.begin(), cmplx.end());
  Rcpp::NumericVector rcppValues(nCmplx);

  typename Rcpp::NumericVector::iterator iValue = rcppValues.begin();
  for (typename Rcpp::List::const_iterator iCmplx = cmplx.begin();
       iCmplx != cmplx.end(); ++iCmplx, ++iValue) {
    Rcpp::NumericVector cmplxVec(*iCmplx);
    typename Rcpp::NumericVector::iterator iCmplxVec = cmplxVec.begin();
    // R is 1-base, while C++ is 0-base
    *iValue = FUNvalues[*iCmplxVec - 1];
    for (; iCmplxVec != cmplxVec.end(); ++iCmplxVec) {
      // R is 1-base, while C++ is 0-base
      *iValue = std::max(*iValue, FUNvalues[*iCmplxVec - 1]);
    }
  }

  // sort
  filtrationSort(rcppCmplx, rcppValues);

  return Rcpp::List::create(rcppCmplx, rcppValues);
}



// RipsFiltration
/** \brief Interface for R code, construct the rips filtration on the input
 *         set of points.
 *
 * @param[out] Rcpp::List     A list
 * @param[in]  X              Either an nxd matrix of coordinates,
 *                            or an nxn matrix of distances of points
 * @param[in]  maxdimension   Max dimension of the homological features to be computed.
 * @param[in]  maxscale       Threshold for the Rips complex
 * @param[in]  dist           "euclidean" for Euclidean distance,
 *                            "arbitrary" for an arbitrary distance
 * @param[in]  library        Either "GUDHI" or "Dionysus"
 * @param[in]  printProgress  Is progress printed?
 * @param[in]  max_num_bars   Write the max_num_pairs most persistent pairs of the
 *                            diagram. Diagram must point to enough memory space for
 *                            3*max_num_pairs double. If there is not enough pairs in the diagram,
 *                            write nothing after.
 */
// [[Rcpp::export]]
Rcpp::List RipsFiltration(
    const Rcpp::NumericMatrix & X,
    const int                   maxdimension,
    const double                maxscale,
    const std::string         & dist,
    const std::string         & library,
    const bool                  printProgress
) {
  if (library[0] == 'G') {
    Gudhi::Simplex_tree<> smplxTree =
      RipsFiltrationGudhi< Gudhi::Simplex_tree<> >(X, maxdimension, maxscale,
        printProgress);
    return filtrationGudhiToRcpp< Rcpp::List, Rcpp::NumericVector >(smplxTree);
  }
  else {

    if (dist[0] == 'e') {
      // RipsDiag for L2 distance
      return filtrationDionysusToRcpp< Rcpp::List, Rcpp::NumericVector >(
        RipsFiltrationDionysus< PairDistances, Generator, FltrR >(X, false,
          maxdimension, maxscale, printProgress));
    }
    else {
      // RipsDiag for arbitrary distance
      return filtrationDionysusToRcpp< Rcpp::List, Rcpp::NumericVector >(
        RipsFiltrationDionysus< PairDistancesA, GeneratorA, FltrRA >(X, true,
          maxdimension, maxscale, printProgress));
    }
  }
}



// RipsDiag
/** \brief Interface for R code, construct the persistence diagram
  * of the Rips complex constructed on the input set of points.
  *
  * @param[out] Rcpp::List     A list
  * @param[in]  X              Either an nxd matrix of coordinates,
  *                            or an nxn matrix of distances of points
  * @param[in]  maxdimension   Max dimension of the homological features to be computed.
  * @param[in]  maxscale       Threshold for the Rips complex
  * @param[in]  dist           "euclidean" for Euclidean distance,
  *                            "arbitrary" for an arbitrary distance
  * @param[in]  libraryFiltration  Either "GUDHI" or "Dionysus"
  * @param[in]  libraryDiag        Either "GUDHI", "Dionysus", or "PHAT"
  * @param[in]  location       Are location of birth point, death point,
  *                            and representative cycles returned?
  * @param[in]  printProgress  Is progress printed?
  * @param[in]  max_num_bars   Write the max_num_pairs most persistent pairs of the
  *                            diagram. Diagram must point to enough memory space for
  *                            3*max_num_pairs double. If there is not enough pairs in the diagram,
  *                            write nothing after.
  */
// [[Rcpp::export]]
Rcpp::List
RipsDiag(const Rcpp::NumericMatrix & X
       , const int                   maxdimension
       , const double                maxscale
       , const std::string         & dist
       , const std::string         & libraryFiltration
       , const std::string         & libraryDiag
       , const bool                  location
       , const bool                  printProgress
) {

	std::vector< std::vector< std::vector< double > > > persDgm;
	std::vector< std::vector< std::vector< unsigned > > > persLoc;
	std::vector< std::vector< std::vector< std::vector< unsigned > > > > persCycle;

	if (libraryFiltration[0] == 'G') {
    Gudhi::Simplex_tree<> st = RipsFiltrationGudhi< Gudhi::Simplex_tree<> >(
        X, maxdimension, maxscale, printProgress);

    // Compute the persistence diagram of the complex
    if (libraryDiag[0] == 'G') {
      int p = 2; //characteristic of the coefficient field for homology
      double min_persistence = 0; //minimal length for persistent intervals
      FiltrationDiagGudhi(
          st, p, min_persistence, maxdimension, printProgress, persDgm);
    }
    else if (libraryDiag[0] == 'D') {
      FltrR filtration = filtrationGudhiToDionysus< FltrR >(st);
      FiltrationDiagDionysus(
          filtration, maxdimension, location, printProgress, persDgm, persLoc,
          persCycle);
    }
    else {
      std::vector< phat::column > cmplx;
      std::vector< double > values;
      phat::boundary_matrix< phat::vector_vector > boundary_matrix;
      filtrationGudhiToPhat< phat::column, phat::dimension >(
        st, cmplx, values, boundary_matrix);
      FiltrationDiagPhat(
        cmplx, values, boundary_matrix, maxdimension, location,
        printProgress, persDgm, persLoc, persCycle);
    }
  }
  else {
    if (dist[0] == 'e') {
      // RipsDiag for L2 distance
      FltrR filtration =
          RipsFiltrationDionysus< PairDistances, Generator, FltrR >(
              X, false, maxdimension, maxscale, printProgress);

      if (libraryDiag[0] == 'D') {
        FiltrationDiagDionysus(
            filtration, maxdimension, location, printProgress, persDgm,
            persLoc, persCycle);
      }
      else if (libraryDiag[0] == 'G') {
        Gudhi::Simplex_tree<> st =
            filtrationDionysusToGudhi< Gudhi::Simplex_tree<> >(filtration);
        int p = 2; //characteristic of the coefficient field for homology
        double min_persistence = 0; //minimal length for persistent intervals
        FiltrationDiagGudhi(
            st, p, min_persistence, maxdimension, printProgress, persDgm);
      }
      else {
        std::vector< phat::column > cmplx;
        std::vector< double > values;
        phat::boundary_matrix< phat::vector_vector > boundary_matrix;
        filtrationDionysusToPhat< phat::column, phat::dimension >(
            filtration, cmplx, values, boundary_matrix);
        FiltrationDiagPhat(
            cmplx, values, boundary_matrix, maxdimension, location,
            printProgress, persDgm, persLoc, persCycle);
      }
    }
    else {
      // RipsDiag for arbitrary distance
      FltrRA filtration =
          RipsFiltrationDionysus< PairDistancesA, GeneratorA, FltrRA >(
              X, true, maxdimension, maxscale, printProgress);

      if (libraryDiag[0] == 'D') {
        FiltrationDiagDionysus(
            filtration, maxdimension, location, printProgress, persDgm,
            persLoc, persCycle);
      }
      else if (libraryDiag[0] == 'G') {
        Gudhi::Simplex_tree<> st =
            filtrationDionysusToGudhi< Gudhi::Simplex_tree<> >(filtration);
        int p = 2; //characteristic of the coefficient field for homology
        double min_persistence = 0; //minimal length for persistent intervals
        FiltrationDiagGudhi(
            st, p, min_persistence, maxdimension, printProgress, persDgm);
      }
      else {
        std::vector< phat::column > cmplx;
        std::vector< double > values;
        phat::boundary_matrix< phat::vector_vector > boundary_matrix;
        filtrationDionysusToPhat< phat::column, phat::dimension >(
          filtration, cmplx, values, boundary_matrix);
        FiltrationDiagPhat(
            cmplx, values, boundary_matrix, maxdimension, location,
            printProgress, persDgm, persLoc, persCycle);
      }
    }
  }

	// Output persistent diagram
	return Rcpp::List::create(
		concatStlToRcpp< Rcpp::NumericMatrix >(persDgm, true, 3),
		concatStlToRcpp< Rcpp::NumericMatrix >(persLoc, false, 2),
		StlToRcppMatrixList< Rcpp::List, Rcpp::NumericMatrix >(persCycle));
}



// AlphaShapeFiltration in GUDHI
/** \brief Interface for R code, construct the persistence diagram of the alpha
 *         shape complex constructed on the input set of points.
 *
 * @param[out] Rcpp::List     A list
 * @param[in]  X              An nx3 matrix of coordinates,
 * @param[in]  printProgress  Is progress printed?
 */
// [[Rcpp::export]]
Rcpp::List AlphaShapeFiltration(
  const Rcpp::NumericMatrix & X,          //points to some memory space
  const bool                  printProgress
) {

  Gudhi::Simplex_tree<> smplxTree =
      AlphaShapeFiltrationGudhi< Gudhi::Simplex_tree<> >(X, printProgress);
  return filtrationGudhiToRcpp< Rcpp::List, Rcpp::NumericVector >(smplxTree);
}



// AlphaShapeDiag in GUDHI
/** \brief Interface for R code, construct the persistence diagram of the alpha
  *        shape complex constructed on the input set of points.
  *
  * @param[out] Rcpp::List     A list
  * @param[in]  X              An nx3 matrix of coordinates,
  * @param[in]  printProgress  Is progress printed?
  */
// [[Rcpp::export]]
Rcpp::List AlphaShapeDiagGudhi(
    const Rcpp::NumericMatrix & X,             //points to some memory space
    const bool                  printProgress
	) {
  std::vector< std::vector< std::vector< double > > > persDgm;
  std::vector< std::vector< std::vector< unsigned > > > persLoc;
  std::vector< std::vector< std::vector< std::vector< unsigned > > > >
      persCycle;

  int coeff_field_characteristic = 2;

  float min_persistence = 0.0;

  Gudhi::Simplex_tree<> simplex_tree =
      AlphaShapeFiltrationGudhi< Gudhi::Simplex_tree<> >(X, printProgress);

  // Compute the persistence diagram of the complex
  FiltrationDiagGudhi(
      simplex_tree, coeff_field_characteristic, min_persistence, 2,
      printProgress, persDgm);

  // Output persistent diagram
  return Rcpp::List::create(
      concatStlToRcpp< Rcpp::NumericMatrix >(persDgm, true, 3),
      concatStlToRcpp< Rcpp::NumericMatrix >(persLoc, false, 2),
      StlToRcppMatrixList< Rcpp::List, Rcpp::NumericMatrix >(persCycle));
}



// AlphaComplexFiltration
/** \brief Interface for R code, construct the persistence diagram of the alpha
 *         complex constructed on the input set of points.
 *
 * @param[out] Rcpp::List     A list
 * @param[in]  X              An nx3 matrix of coordinates,
 * @param[in]  maxalphasquare Threshold for the Alpha complex,
 * @param[in]  printProgress  Is progress printed?
 */
// [[Rcpp::export]]
Rcpp::List AlphaComplexFiltration(
  const Rcpp::NumericMatrix & X,             //points to some memory space
  const bool                  printProgress
) {

  Gudhi::Simplex_tree<> smplxTree =
      AlphaComplexFiltrationGudhi< Gudhi::Simplex_tree<> >(X, printProgress);
  return filtrationGudhiToRcpp< Rcpp::List, Rcpp::NumericVector >(smplxTree);
}



// AlphaComplexDiag in GUDHI
/** \brief Interface for R code, construct the persistence diagram of the alpha
  *        complex constructed on the input set of points.
  *
  * @param[out] Rcpp::List     A list
  * @param[in]  X              An nx3 matrix of coordinates,
  * @param[in]  maxalphasquare Threshold for the Alpha complex,
  * @param[in]  printProgress  Is progress printed?
  */
// [[Rcpp::export]]
Rcpp::List AlphaComplexDiagGudhi(
    const Rcpp::NumericMatrix & X,             //points to some memory space
    const bool                  printProgress
	) {
  std::vector< std::vector< std::vector< double > > > persDgm;
  std::vector< std::vector< std::vector< unsigned > > > persLoc;
  std::vector< std::vector< std::vector< std::vector< unsigned > > > >
      persCycle;

  using Kernel = CGAL::Epick_d< CGAL::Dynamic_dimension_tag>;
  using Point = Kernel::Point_d;

  int coeff_field_characteristic = 2;

  float min_persistence = 0.0;

  Gudhi::Simplex_tree<> alpha_complex_from_points =
      AlphaComplexFiltrationGudhi< Gudhi::Simplex_tree<> >(X, printProgress);

  // Compute the persistence diagram of the complex
  FiltrationDiagGudhi(
      alpha_complex_from_points, coeff_field_characteristic, min_persistence,
      2, printProgress, persDgm);

  // Output persistent diagram
  return Rcpp::List::create(
      concatStlToRcpp< Rcpp::NumericMatrix >(persDgm, true, 3),
      concatStlToRcpp< Rcpp::NumericMatrix >(persLoc, false, 2),
      StlToRcppMatrixList< Rcpp::List, Rcpp::NumericMatrix >(persCycle));
}