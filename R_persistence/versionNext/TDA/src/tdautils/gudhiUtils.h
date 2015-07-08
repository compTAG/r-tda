#include <R.h>
#include <R_ext/Print.h>

#include <gudhi/reader_utils.h>
#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Persistent_cohomology/Field_Zp.h>

typedef int        Vertex_handle;
typedef double     Filtration_value;



template<typename SimplexTree>
void computePersistenceGUDHI(SimplexTree& simplexTree,
		const unsigned coeffFieldCharacteristic, const double minPersistence,
		const unsigned maxdimension,
		std::vector< std::vector< std::vector< double > > > &persDgm) {
	// Compute the persistence diagram of the complex
	Gudhi::persistent_cohomology::Persistent_cohomology< Gudhi::Simplex_tree<>, Gudhi::persistent_cohomology::Field_Zp > pcoh(simplexTree);
	pcoh.init_coefficients(coeffFieldCharacteristic); //initializes the coefficient field for homology

	pcoh.compute_persistent_cohomology(minPersistence); //compute persistent homology

	
	typename Gudhi::persistent_cohomology::Persistent_cohomology< Gudhi::Simplex_tree<>, Gudhi::persistent_cohomology::Field_Zp >::cmp_intervals_by_length cmp(pcoh.cpx_);
	pcoh.persistent_pairs_.sort(cmp);
	persDgm.resize(maxdimension + 1);
	std::vector< double > persDgmPoint(2);
	for (auto pair : pcoh.persistent_pairs_) {
		persDgmPoint[0] = pcoh.cpx_->filtration(get<0>(pair));
		persDgmPoint[1] = pcoh.cpx_->filtration(get<1>(pair));
		persDgm[pcoh.cpx_->dimension(get<0>(pair))].push_back(persDgmPoint);
	}

	// write diagram on the output file
	//	pcoh.write_output_diagram(diagram_name);

	// or write the most persistent points in diagram (max_num_bars should be an input) 
	//  pcoh.most_persistent_bars(diagram, *max_num_bars);
}
