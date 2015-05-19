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

// for GUDHI
#include <tdautils/gudhiUtils.h>

// for Dionysus
#include <tdautils/dionysusUtils.h>

// for phat
#include <tdautils/phatUtils.h>



/**
 * grid function by Brittany T. Fasy
 * modified by Jisu Kim for
 * arbitrary dimension & using memory as an input & setting maximum dimension.
 */
// [[Rcpp::export]]
Rcpp::List GridDiag(const Rcpp::NumericVector& FUNvalues,
		const Rcpp::IntegerVector& gridDim, const int maxdimension,
		const std::string& decomposition, const std::string& library,
		const bool location, const bool printProgress) {
#ifdef LOGGING
	//rlog::RLogInit(argc, argv);

	stdoutLog.subscribeTo(RLOG_CHANNEL("topology/persistence"));
	//stdoutLog.subscribeTo(RLOG_CHANNEL("topology/chain"));
	//stdoutLog.subscribeTo(RLOG_CHANNEL("topology/vineyard"));
#endif

	Fltr f;

	// Generate simplicial complex from function values and grid
	if (decomposition[0] == '5') {
		simplicesFromGrid(f, FUNvalues, gridDim, maxdimension + 1);
	}
	if (decomposition[0] == 'b') {
		simplicesFromGridBarycenter(f, FUNvalues, gridDim, maxdimension + 1);
	}
	if (printProgress) {
		Rprintf("# Generated complex of size: %d \n", f.size());
	}

	// Sort simplicial complex
	f.sort(Smplx::DataComparison());

	// Compute persistent homology from sorted simplicial complex
	std::vector< std::vector< std::vector< double > > > persDgm;
	std::vector< std::vector< std::vector< unsigned > > > persLoc;
	std::vector< std::vector< std::set< unsigned > > > persCycle;

	if (library[0] == 'D') {
		persistentPairsDionysus(f, maxdimension, FUNvalues, location,
				printProgress, persDgm, persLoc, persCycle);
	}
	if (library[0] == 'P') {
		 persistentPairsPhat(f, maxdimension, FUNvalues, location,
			 printProgress, persDgm, persLoc);
	}

	// Output persistent diagram
	return Rcpp::List::create(
		concatStlToRcpp< Rcpp::NumericMatrix >(persDgm, true, 3),
		concatStlToRcpp< Rcpp::NumericMatrix >(persLoc, false, 2),
		StlToRcppList< Rcpp::List, Rcpp::NumericVector >(persCycle));


	//// Output persistent diagram
	//std::ofstream outfile;
	//outfile.open("outputTDA.txt");
	//for (int dgmsIdx = 0; dgmsIdx <= maxdimension; ++dgmsIdx)
	//{
	//	outfile << dgmsIdx << std::endl;
	//	std::vector< std::vector< double > >::const_iterator persdgmIdx;
	//	for (persdgmIdx = persDgm[dgmsIdx].begin(); persdgmIdx != persDgm[dgmsIdx].end(); ++persdgmIdx) {
	//		outfile << (*persdgmIdx)[0] << " " << (*persdgmIdx)[1] << std::endl;
	//	}
	//	outfile << std::endl;
	//}
	//if (location)
	//{
	//	outfile << "Location" << std::endl;
	//	for (int dgmsIdx = 0; dgmsIdx <= maxdimension; ++dgmsIdx) {
	//		std::vector< std::vector< unsigned int > >::const_iterator persLocIdx;
	//		for (persLocIdx = persLoc[dgmsIdx].begin(); persLocIdx != persLoc[dgmsIdx].end(); ++persLocIdx) {
	//			outfile << (*persLocIdx)[0] << " " << (*persLocIdx)[1] << std::endl;
	//		}
	//	}

	//	outfile << std::endl << "Cycle" << std::endl;
	//	for (int dgmsIdx = 0; dgmsIdx <= maxdimension; ++dgmsIdx) {
	//		std::vector< std::set< unsigned int > >::const_iterator persCycleIdx;
	//		std::set< unsigned int >::const_iterator persCyclePointIdx;
	//		for (persCycleIdx = persCycle[dgmsIdx].begin(); persCycleIdx != persCycle[dgmsIdx].end(); ++persCycleIdx) {
	//			for (persCyclePointIdx = persCycleIdx->begin(); persCyclePointIdx != persCycleIdx->end(); ++persCyclePointIdx) {
	//				outfile << (*persCyclePointIdx) << " ";
	//			}
	//			outfile << std::endl;
	//		}
	//	}
	//}
	//outfile.close();
}



extern "C" {

	void rips(int* dimInput, double* maxInput, int* printInput)
	{
		bool printstatus=printInput[0];
		Dimension               skeleton;
		DistanceType            max_distance;
		std::string             infilename, diagram_name;
		
		infilename = "inputTDA.txt";
		diagram_name="outputTDA.txt";
		skeleton= dimInput[0];
		max_distance= maxInput[0];
		
		std::ofstream           diagram_out(diagram_name.c_str());
		
		/*if (printstatus){
			Rprintf("Diagram %s \n", diagram_name.c_str());
			//std::cout << "Diagram:         " << diagram_name << std::endl;
		}*/
		PointContainer          points;
		read_points(infilename, points);

		PairDistances           distances(points);
		Generator               rips(distances);
		Generator::Evaluator    size(distances);
		FltrR                    f;
	
		// Generate 2-skeleton of the Rips complex for epsilon = 50
		rips.generate(skeleton, max_distance, make_push_back_functor(f));
		if (printstatus){
			Rprintf("# Generated complex of size: %d \n", f.size());
			//std::cout << "# Generated complex of size: " << f.size() << std::endl;
		}
		// Generate filtration with respect to distance and compute its persistence
		f.sort(Generator::Comparison(distances));

		Timer persistence_timer; persistence_timer.start();
		PersistenceR p(f);
		p.pair_simplices(printstatus);
		persistence_timer.stop();

	#if 1
		// Output cycles
		PersistenceR::SimplexMap<FltrR>   m = p.make_simplex_map(f);
		for (PersistenceR::iterator cur = p.begin(); cur != p.end(); ++cur)
		{
//			const PersistenceR::Cycle& cycle = cur->cycle;

			if (!cur->sign())        // only negative simplices have non-empty cycles
			{
				PersistenceR::OrderIndex birth = cur->pair;      // the cycle that cur killed was born when we added birth (another simplex)

				const SmplxR& b = m[birth];
				const SmplxR& d = m[cur];
			
				// if (b.dimension() != 1) continue;
				// std::cout << "Pair: (" << size(b) << ", " << size(d) << ")" << std::endl;
				if (b.dimension() >= skeleton) continue;
				diagram_out << b.dimension() << " " << size(b) << " " << size(d) << std::endl;
			} else if (cur->unpaired())    // positive could be unpaired
			{
				const SmplxR& b = m[cur];
				// if (b.dimension() != 1) continue;
			
				// std::cout << "Unpaired birth: " << size(b) << std::endl;
				// cycle = cur->chain;      // TODO
				if (b.dimension() >= skeleton) continue;
				diagram_out << b.dimension() << " " << size(b) << " inf" << std::endl;
			}

			// Iterate over the cycle
			// for (PersistenceR::Cycle::const_iterator si =  cycle.begin();
			//                                                          si != cycle.end();     ++si)
			// {
			//     const SmplxR& s = m[*si];
			//     //std::cout << s.dimension() << std::endl;
			//     const SmplxR::VertexContainer& vertices = s.vertices();          // std::vector<Vertex> where Vertex = Distances::IndexType
			//     AssertMsg(vertices.size() == s.dimension() + 1, "dimension of a simplex is one less than the number of its vertices");
			//     std::cout << vertices[0] << " " << vertices[1] << std::endl;
			// }
		}
	#endif
	
		if (printstatus){	
			persistence_timer.check("# Persistence timer");
		}

	}



	void ripsArbit(int* dimInput, double* maxInput, int* printInput)
	{
		bool printstatus=printInput[0];
		Dimension               skeleton;
		DistanceTypeA            max_distance;
		std::string             infilename, diagram_name;

		infilename = "inputTDA.txt";
		diagram_name="outputTDA.txt";
		skeleton= dimInput[0];
		max_distance= maxInput[0];
		
		std::ofstream           diagram_out(diagram_name.c_str());
		/*if (printstatus){
			Rprintf("Diagram %s \n", diagram_name.c_str());			
			//std::cout << "Diagram:         " << diagram_name << std::endl;
		}*/
		PointContainer          points;
		read_points2(infilename, points);

		PairDistancesA           distances(points);
		GeneratorA               rips(distances);
		GeneratorA::Evaluator    size(distances);
		FltrRA                    f;
	
		// Generate 2-skeleton of the Rips complex for epsilon = 50
		rips.generate(skeleton, max_distance, make_push_back_functor(f));
		if (printstatus){
			Rprintf("# Generated complex of size: %d \n", f.size());
			//std::cout << "# Generated complex of size: " << f.size() << std::endl;
		}
		// Generate filtration with respect to distance and compute its persistence
		f.sort(GeneratorA::Comparison(distances));

		Timer persistence_timer; persistence_timer.start();
		PersistenceR p(f);
		p.pair_simplices(printstatus);
		persistence_timer.stop();

	#if 1
		// Output cycles
		PersistenceR::SimplexMap<FltrRA>   m = p.make_simplex_map(f);
		for (PersistenceR::iterator cur = p.begin(); cur != p.end(); ++cur)
		{
//			const PersistenceR::Cycle& cycle = cur->cycle;

			if (!cur->sign())        // only negative simplices have non-empty cycles
			{
				PersistenceR::OrderIndex birth = cur->pair;      // the cycle that cur killed was born when we added birth (another simplex)

				const SmplxRA& b = m[birth];
				const SmplxRA& d = m[cur];
			
				// if (b.dimension() != 1) continue;
				// std::cout << "Pair: (" << size(b) << ", " << size(d) << ")" << std::endl;
				if (b.dimension() >= skeleton) continue;
				diagram_out << b.dimension() << " " << size(b) << " " << size(d) << std::endl;
			} else if (cur->unpaired())    // positive could be unpaired
			{
				const SmplxRA& b = m[cur];
				// if (b.dimension() != 1) continue;
			
				// std::cout << "Unpaired birth: " << size(b) << std::endl;
				// cycle = cur->chain;      // TODO
				if (b.dimension() >= skeleton) continue;
				diagram_out << b.dimension() << " " << size(b) << " inf" << std::endl;
			}

			// Iterate over the cycle
			// for (PersistenceR::Cycle::const_iterator si =  cycle.begin();
			//                                                          si != cycle.end();     ++si)
			// {
			//     const SmplxRA& s = m[*si];
			//     //std::cout << s.dimension() << std::endl;
			//     const SmplxR::VertexContainer& vertices = s.vertices();          // std::vector<Vertex> where Vertex = Distances::IndexType
			//     AssertMsg(vertices.size() == s.dimension() + 1, "dimension of a simplex is one less than the number of its vertices");
			//     std::cout << vertices[0] << " " << vertices[1] << std::endl;
			// }
		}
	#endif
		if (printstatus){	
			persistence_timer.check("# Persistence timer");
		}
	}

}



// [[Rcpp::export]]
double Bottleneck(
		const Rcpp::NumericMatrix& Diag1, const Rcpp::NumericMatrix& Diag2) {
	return bottleneck_distance(RcppToDionysus< PersistenceDiagram<> >(Diag1),
			RcppToDionysus< PersistenceDiagram<> >(Diag2));
}



// [[Rcpp::export]]
double Wasserstein(const Rcpp::NumericMatrix& Diag1,
		const Rcpp::NumericMatrix& Diag2, const int p) {
	return wasserstein_distance(RcppToDionysus< PersistenceDiagram<> >(Diag1),
			RcppToDionysus< PersistenceDiagram<> >(Diag2), p);
}



// KDE function on a Grid
// [[Rcpp::export]]
Rcpp::NumericVector Kde(const Rcpp::NumericMatrix& X,
		const Rcpp::NumericMatrix& Grid, const double h,
		const Rcpp::NumericVector& weight, const bool printProgress) {
	const double pi = 3.141592653589793;
	const unsigned dimension = Grid.ncol();
	const unsigned gridNum = Grid.nrow();
	const double den = pow(h, (int)dimension) * pow(2 * pi, dimension / 2.0);
	Rcpp::NumericVector kdeValue(gridNum);

	int counter = 0;
	const int totalCount = gridNum;
	int percentageFloor = 0;

	if (printProgress) {
		printProgressFrame(Rprintf);

		for (unsigned gridIdx = 0; gridIdx < gridNum; ++gridIdx) {
			kdeValue[gridIdx] = oneKernel(matrixRow< std::vector< double > >(
					Grid, gridIdx), X, h, weight) / den;
						
			//printProgress
			printProgressAmount(Rprintf, counter, totalCount, percentageFloor);
		}
		Rprintf("\n");

	} else { //no printProgress
		for (unsigned gridIdx = 0; gridIdx < gridNum; ++gridIdx) {
			kdeValue[gridIdx] = oneKernel(matrixRow< std::vector< double > >(
					Grid, gridIdx), X, h, weight) / den;
		}
	}

	return (kdeValue);
}



// kernel Dist function on a Grid
// [[Rcpp::export]]
Rcpp::NumericVector KdeDist(const Rcpp::NumericMatrix& X,
		const Rcpp::NumericMatrix& Grid, const double h,
		const Rcpp::NumericVector& weight, const bool printProgress) {
	const unsigned sampleNum = X.nrow();
	const unsigned gridNum = Grid.nrow();
	// first = sum K_h(X_i, X_j), second = K_h(x, x), third = sum K_h(x, X_i)
	double first = 0.0;
	const double second = 1.0;
	double third;
	Rcpp::NumericVector kdeDistValue(gridNum);

	int counter = 0;
	const int totalCount = sampleNum + gridNum;
	int percentageFloor = 0;

	if (printProgress)
	{
		printProgressFrame(Rprintf);

		for (unsigned sampleIdx = 0; sampleIdx < sampleNum; ++sampleIdx) {
			first += oneKernel(
					matrixRow< std::vector< double > >(X, sampleIdx), X, h, weight);

			// printProgress
			printProgressAmount(Rprintf, counter, totalCount, percentageFloor);
		}
		first /= sampleNum;

		for (unsigned gridIdx = 0; gridIdx < gridNum; ++gridIdx) {
			third = oneKernel(
					matrixRow< std::vector< double > >(Grid, gridIdx), X, h, weight);
			kdeDistValue[gridIdx] = std::sqrt(first + second - 2 * third);

			// printProgress
			printProgressAmount(Rprintf, counter, totalCount, percentageFloor);
		}
		Rprintf("\n");		

	} else { //no printProgress
		for (unsigned sampleIdx = 0; sampleIdx < sampleNum; ++sampleIdx) {
			first += oneKernel(
					matrixRow< std::vector< double > >(X, sampleIdx), X, h, weight);
		}
		first /= sampleNum;

		for (unsigned gridIdx = 0; gridIdx < gridNum; ++gridIdx) {
			third = oneKernel(
					matrixRow< std::vector< double > >(Grid, gridIdx), X, h, weight);
			// first = sum K_h(X_i, X_j), second = K_h(x, x), third = sum K_h(x, X_i)
			kdeDistValue[gridIdx] = std::sqrt(first + second - 2 * third);
		}
	}

	return (kdeDistValue);
}



// distance to measure function on a Grid
// [[Rcpp::export]]
Rcpp::NumericVector Dtm(const Rcpp::NumericMatrix& knnIndex,
		const Rcpp::NumericMatrix& knnDistance,
		const Rcpp::NumericVector& weight, const double weightBound) {
	const unsigned gridNum = knnIndex.nrow();
	double distanceTemp, weightTemp, weightSumTemp;
	Rcpp::NumericVector dtmValue(gridNum, 0.0);

	for (unsigned gridIdx = 0; gridIdx < gridNum; ++gridIdx) {
		weightSumTemp = 0.0;
		for (unsigned kIdx = 0; weightSumTemp < weightBound; ++kIdx) {
			distanceTemp = knnDistance[gridIdx + kIdx * gridNum];
			weightTemp = std::min(weight[knnIndex[gridIdx + kIdx * gridNum] - 1],
					weightBound - weightSumTemp);
			weightSumTemp += weightTemp;
			dtmValue[gridIdx] += distanceTemp * distanceTemp * weightTemp;
		}
		dtmValue[gridIdx] = std::sqrt(dtmValue[gridIdx] / weightBound);
	}
	return (dtmValue);
}


	extern "C" {

	// GUDHI RIPS
	/** \brief Interface for R code, construct the persistence diagram 
	  * of the Rips complex constructed on the input set of points.
	  *
	  * @param[out] void            every function called by R must return void
	  * @param[in]  point           pointer toward the coordinates of all points. Format 
	  *                             must be X11 X12 ... X1d X21 X22 ... X2d X31 ...
	  * @param[in]  dim             embedding dimension 
	  * @param[in]  num_points      number of points. The input point * must be a 
	  *                             pointer toward num_points*dim double exactly.
	  * @param[in]  rips_threshold  threshold for the Rips complex
	  * @param[in]  max_complex_dim maximal dimension of the Rips complex
	  * @param[in]  diagram         where to output the diagram. The format must be dimension birth death.
	  * @param[in]  max_num_bars    write the max_num_pairs most persistent pairs of the
	  *                             diagram. diagram must point to enough memory space for
	  *                             3*max_num_pairs double. If there is not enough pairs in the diagram,
	  *                             write nothing after.
	  */
	void 
	rips_persistence_diagram_GUDHI( double * points          //points to some memory space
							, int    * dim
							, int    * num_points
							, double * rips_threshold
							, int    * max_complex_dim
							, double * diagram         //points to some memory space
							, int    * printInput)

	{
	  bool printstatus=printInput[0];  
	  std::string diagram_name="outputTDA.txt";


	  // Turn the input points into a range of points
	  typedef std::vector<double> Point_t;
	  std::vector< Point_t > point_set(*num_points, Point_t(*dim));
	  double * curr_coordinate = points;
	  for(int i = 0; i < *num_points; ++i) {
		for(int j = 0; j < *dim; ++j) {
		  point_set[i][j] = *curr_coordinate;
		  ++curr_coordinate;
		}
	  }
	// Compute the proximity graph of the points
	  Graph_t prox_graph = compute_proximity_graph( point_set, *rips_threshold
												  , euclidean_distance<Point_t> );

	// Construct the Rips complex in a Simplex Tree
	// Construct the Rips complex in a Simplex Tree
	  Simplex_tree<> st;        
	  st.insert_graph(prox_graph); // insert the proximity graph in the simplex tree
	  st.expansion( *max_complex_dim ); // expand the graph until dimension dim_max

		if (printstatus){
			Rprintf("# Generated complex of size: %d \n", st.num_simplices());
			// std::cout << st.num_simplices() << " simplices \n";
		}

	// Sort the simplices in the order of the filtration
	  st.initialize_filtration();

	// Compute the persistence diagram of the complex
	  int p = 2; //characteristic of the coefficient field for homology
	  double min_persistence = 0; //minimal length for persistent intervals
	  Persistent_cohomology< Simplex_tree<>, Field_Zp > pcoh( st );
	  pcoh.init_coefficients( p ); //initilizes the coefficient field for homology
	  pcoh.compute_persistent_cohomology( min_persistence ); //compute persistent homology

	// write diagram on the output file
	  pcoh.write_output_diagram(diagram_name);
	
	// or write the most persistent points in diagram (max_num_bars should be an input) 
	//  pcoh.most_persistent_bars(diagram, *max_num_bars);
	}



} //end extern
