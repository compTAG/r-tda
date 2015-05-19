#ifndef __DIONYSUSUTILS_H__
#define __DIONYSUSUTILS_H__

#include <topology/simplex.h>
#include <string>
#include <sstream>
#include <cstdlib>



std::vector<unsigned int> getVertices(const Simplex<unsigned, double>& smp)
{
	std::stringstream sstr;
	std::vector<unsigned int> vertices(smp.dimension() + 1);
	sstr << smp;
	std::string vtxStr;
	std::getline(sstr,vtxStr,'<');
	unsigned int vtxIdx;
	for (vtxIdx = 0; vtxIdx < (unsigned)smp.dimension(); ++vtxIdx)
	{
		std::getline(sstr, vtxStr, ',');
		vertices[vtxIdx] = (unsigned)std::atoi(vtxStr.c_str());
	}
	std::getline(sstr,vtxStr,'>');
	vertices[vtxIdx] = (unsigned)std::atoi(vtxStr.c_str());
	return vertices;
}



template<typename RealVector>
unsigned int getLocation(
		const Simplex<unsigned, double>& smp, const RealVector& FUNvalues)
{
	std::vector< unsigned > vertices;
	std::vector< unsigned >::const_iterator vertexItr;
	unsigned vertex;
	vertices = getVertices(smp);
	vertex = *(vertices.begin());
	for (vertexItr = vertices.begin(); vertexItr != vertices.end(); ++vertexItr)
	{
		if (FUNvalues[*vertexItr] > FUNvalues[vertex])
		{
			vertex = *vertexItr;
		}
	}
	return vertex + 1;
}



template<typename Flt, typename RealVector>
void persistentPairsDionysus(Flt f, const unsigned maxdimension,
		const RealVector& FUNvalues, const bool location, const bool printProgress,
		std::vector< std::vector< std::vector< double > > > &persDgm,
		std::vector< std::vector< std::vector< unsigned int > > > &persLoc,
		std::vector< std::vector< std::set< unsigned int > > >& persCycle) {
	Timer persistence_timer;
	persistence_timer.start();

	// Compute persistent homology from sorted simplicial complex
	Persistence p(f); // initialize persistence
	if (!location) {
		p.pair_simplices(printProgress); // pair simplices
	}
	else {
		if (printProgress) {
			p.pair_simplices(p.begin(), p.end(), true, Persistence::PairVisitor(p.size()));
		}
		else {
			p.pair_simplices(p.begin(), p.end(), true, Persistence::PairVisitorNoProgress());
		}
	}

	persistence_timer.stop();

	// TODO: why doesn't this work? rLog(rlmain, "testing");  

	Persistence::SimplexMap<Fltr> m = p.make_simplex_map(f);
	std::map<Dimension, PersistenceDiagram<> > dgms;
	init_diagrams(dgms, p.begin(), p.end(),
		evaluate_through_map(m, Smplx::DataEvaluator()),
		evaluate_through_map(m, Smplx::DimensionExtractor()));


	// Save persistent diagram & trace back birth & death simplex
	persDgm.resize(maxdimension + 1);
	std::vector< double > persDgmPoint(2);
	std::map<Dimension, PersistenceDiagram<> >::const_iterator dgmItr;
	persDgm.resize(maxdimension + 1);
	for (dgmItr = dgms.begin(); dgmItr != dgms.end(); ++dgmItr) {
		PersistenceDiagram<>::const_iterator ptItr;
		for (ptItr = (dgmItr->second).begin(); ptItr != (dgmItr->second).end();
			++ptItr) {
			persDgmPoint[0] = ptItr->x();
			persDgmPoint[1] = ptItr->y();
			persDgm[dgmItr->first].push_back(persDgmPoint);
		}
	}

	if (location) {
		persLoc.resize(maxdimension + 1);
		persCycle.resize(maxdimension + 1);
		std::vector< unsigned int > persLocPoint(2);
		std::set< unsigned int > persCyclePoint;
		for (Persistence::iterator cur = p.begin(); cur != p.end(); ++cur) {
			// positive simplices corresponds to
			// negative simplices having non-empty cycles
			if (cur->sign()) {
				// first consider that cycle is paired
				if (!cur->unpaired()) {
					// the cycle that was born at cur is killed 
					// when we added death (another simplex)
					Persistence::OrderIndex death = cur->pair;

					const Smplx& b = m[cur];
					const Smplx& d = m[death];
					if (b.data() < d.data() && b.dimension() <= maxdimension) {
						persLocPoint[0] = getLocation(b, FUNvalues);
						persLocPoint[1] = getLocation(d, FUNvalues);
						persLoc[b.dimension()].push_back(persLocPoint);

						// Iterate over the cycle
						persCyclePoint.clear();
						const Persistence::Cycle& cycle = death->cycle;
						for (Persistence::Cycle::const_iterator si = cycle.begin();
								si != cycle.end(); ++si) {
							const Smplx& s = m[*si];
							const Smplx::VertexContainer& vertices = s.vertices();    // std::vector<Vertex> where Vertex = Distances::IndexType
							Smplx::VertexContainer::const_iterator vtxItr;
							for (vtxItr = vertices.begin(); vtxItr != vertices.end();
									++vtxItr) {
								persCyclePoint.insert(*vtxItr + 1);
							}
						}
						persCycle[b.dimension()].push_back(persCyclePoint);
					}
				}
				else {    // cycles can be unpaired
					const Smplx& b = m[cur];
					persLocPoint[0] = getLocation(b, FUNvalues);
					persLocPoint[1] = (unsigned)(std::max_element(
							FUNvalues.begin(), FUNvalues.end()) - FUNvalues.begin() + 1);
					persLoc[b.dimension()].push_back(persLocPoint);

					// Iterate over the cycle
					persCyclePoint.clear();
					persCycle[b.dimension()].push_back(persCyclePoint);
				}
			}
		}
	}

	if (printProgress) {
		persistence_timer.check("# Persistence timer");
	}
}



# endif // __DIONYSUSUTILS_H__