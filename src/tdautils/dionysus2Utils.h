#ifndef __DIONYSUS2UTILS_H__
#define __DIONYSUS2UTILS_H__

#include <dionysus/simplex.h>
#include <dionysus/filtration.h>
#include <dionysus/fields/z2.h>
#include <dionysus/ordinary-persistence.h>
#include <dionysus/standard-reduction.h>
#include <dionysus/row-reduction.h>
#include <dionysus/diagram.h>
#include "diagramDS.h"
#include <utilities/timer.h>

namespace d = dionysus;

// FiltrationDiag in Dionysus2
/** \brief Construct the persistence diagram from the filtration using library
*         Dionysus.
*
* @param[out] void           Void
* @param[in]  filtration     The input filtration
* @param[in]  maxdimension   Max dimension of the homological features to be
*                            computed
* @param[in]  location       Are location of birth point, death point, and
*                            representative cycles returned?
* @param[in]  printProgress  Is progress printed?
* @param[in]  persDgm        Memory space for the resulting persistence
*                            diagram
* @param[in]  persLoc        Memory space for the resulting birth points and
*                            death points
* @param[in]  persCycle      Memory space for the resulting representative
*                            cycles
* @param[in]  max_num_bars   Write the max_num_pairs most persistent pairs of the
*                            diagram. Diagram must point to enough memory space for
*                            3*max_num_pairs double. If there is not enough pairs in the diagram,
*                            write nothing after.
*/


/**
 * Class: EvaluatePushBack<Container>
 *
 * Push back the simplex and the evaluated value
 */
template< typename Container, typename Evaluator >
class EvaluatePushBack2 {

public:
  EvaluatePushBack2(Container & argContainer, const Evaluator & argEvaluator) :
    container(argContainer), evaluator(argEvaluator) {}

  void operator()(const typename Container::value_type & argSmp) const {
    typename Container::value_type smp(argSmp.dimension(),argSmp.begin(),argSmp.end(), evaluator(argSmp));
    container.push_back(smp);
  }

private:
  Container & container;
  const Evaluator & evaluator;
};

template< typename VertexList, typename Evaluator >
unsigned getLocation(const VertexList & vertices, const Evaluator & evaluator) {
  typename VertexList::const_iterator vertexItr;
  unsigned vertex = *(vertices.begin());
	for (vertexItr = vertices.begin(); vertexItr != vertices.end(); ++vertexItr) {
		if (evaluator[*vertexItr] > evaluator[vertex]) {
			vertex = *vertexItr;
		}
	}
	return vertex + 1;
}

template< typename Simplex, typename Locations, typename Cycles,
          typename Persistence, typename Evaluator, typename SimplexMap,
          typename Filtration >
inline void initLocations(
    Locations & locations, Cycles & cycles, const Persistence & p,
    const Evaluator & evaluator, const SimplexMap & m,
    const unsigned maxdimension, const Filtration & filtration) {

	unsigned verticesMax = 0;
	for (typename Filtration::OrderConstIterator iFltr = filtration.begin();
       iFltr != filtration.end(); ++iFltr) {
		const typename Filtration::Simplex & c = *(iFltr);
		if (c.dimension() == 0) {
			verticesMax = std::max(verticesMax, *(c.begin()));
		}
	}

  // vertices range from 0 to verticesMax
  std::vector< double > verticesValues(
      verticesMax + 1, -std::numeric_limits< double >::infinity());

  for (typename Filtration::OrderConstIterator iFltr = filtration.begin();
       iFltr != filtration.end(); ++iFltr) {
		const typename Filtration::Cell & c = *(iFltr);
		if(c.dimension() == 0) {
			verticesValues[*(c.begin())] = c.data();
		}
	}

	locations.resize(maxdimension + 1);
	cycles.resize(maxdimension + 1);
	typename Locations::value_type::value_type persLocPoint(2);
	typename Cycles::value_type::value_type persBdy;
	typename Cycles::value_type::value_type::value_type persSimplex;
	for (typename Persistence::iterator cur = p.begin(); cur != p.end(); ++cur) {
		// positive simplices corresponds to
		// negative simplices having non-empty cycles
		if (cur->sign()) {
			// first consider that cycle is paired
			if (!cur->unpaired()) {
				// the cycle that was born at cur is killed 
				// when we added death (another simplex)
				const typename SimplexMap::key_type& death = cur->pair;

				//const typename SimplexMap::value_type& b = m[cur];
				//const typename SimplexMap::value_type& d = m[death];
        const typename Filtration::Cell & b = m[cur];
        const typename Filtration::Cell & d = m[death];
				if ((unsigned)b.dimension() > maxdimension) {
					continue;
				}
				if (evaluator(b) < evaluator(d)) {
					persLocPoint[0] = getLocation(b.vertices(), verticesValues);
					persLocPoint[1] = getLocation(d.vertices(), verticesValues);
					locations[b.dimension()].push_back(persLocPoint);

					// Iterate over the cycle
					persBdy.clear();
					const typename Persistence::Cycle& cycle = death->cycle;
					for (typename Persistence::Cycle::const_iterator
						si = cycle.begin(); si != cycle.end(); ++si) {
						persSimplex.clear();
						const typename Simplex::VertexContainer&
							vertices = m[*si].vertices();    // std::vector<Vertex> where Vertex = Distances::IndexType
						typename Simplex::VertexContainer::const_iterator vtxItr;
						for (vtxItr = vertices.begin(); vtxItr != vertices.end();
							++vtxItr) {
							persSimplex.push_back(*vtxItr + 1);
						}
						persBdy.push_back(persSimplex);
					}
					cycles[b.dimension()].push_back(persBdy);
				}
			}
			else {    // cycles can be unpaired
				const typename SimplexMap::value_type& b = m[cur];
				if ((unsigned)b.dimension() > maxdimension) {
					continue;
				}
				persLocPoint[0] = getLocation(b.vertices(), verticesValues);
				persLocPoint[1] = (unsigned)(
            std::max_element(verticesValues.begin(), verticesValues.end())
            - verticesValues.begin() + 1);
				locations[b.dimension()].push_back(persLocPoint);

				// Iterate over the cycle
				persBdy.clear();
				cycles[b.dimension()].push_back(persBdy);
			}
		}
	}
}

template <typename Persistence, typename Filtration>
void FiltrationDiagDionysus2(
        const Filtration &filtration,
        const int maxdimension,
        const bool location,
        const bool printProgress,
        std::vector< std::vector< std::vector< double > > >   & persDgm,
        std::vector< std::vector< std::vector< unsigned > > > & persLoc,
        std::vector< std::vector< std::vector< std::vector< unsigned > > > > & persCycle
) {
    
    // Assume that Persistence that is passed in is Persistence2, Filtration is Fltr2
    // Create and Calculate Persistence
    d::Z2Field k;
    d::RowReduction<d::Z2Field> reduce(k);
    reduce(filtration);
    
    typedef decltype(reduce.persistence().pair(0)) Index;
    typedef float Value;
    // Putting persistence into diagrams data structure
    diagramDS::DiagramDS<Value, Index> dgms(
            reduce.persistence(), 
            filtration, 
            [&](const Smplx2& s) -> float  { return filtration.index(s);}, 
            [](typename Persistence::Index i) { return i; }
        );
    // Fill in Diagram
    persDgm.resize(dgms.getDiagrams().size());
    Index _ = 0;
    for (auto &dgm : dgms.getDiagrams())
    {
        for (auto &pt : dgm)
        {
            std::vector<double> pt_;
            if (pt.death() == std::numeric_limits<Value>::infinity()) {
                pt_ = {filtration[pt.birth()].data(),pt.death()};
            } else {
                pt_ = {filtration[pt.birth()].data(),filtration[pt.death()].data()};
            }
            if (pt_[0] != pt_[1])
                persDgm[_].push_back(pt_);
        }
        _++;   
    }
    // Capping at maxdimension
    if (persDgm.size() > maxdimension) {
        persDgm.resize(maxdimension + 1);
    }
    if (location) {
    initLocations< typename Filtration::Cell >(
        persLoc, persCycle, p, typename Filtration::Cell::DataEvaluator(),
        m, maxdimension, filtration);
    }
}
template< typename Distances, typename Generator, typename Filtration,
          typename RealMatrix, typename Print >
inline Filtration RipsFiltrationDionysus2(
    const RealMatrix & X,
    const unsigned     nSample, 
    const unsigned     nDim,
    const bool         is_row_names,
    const int          maxdimension,
    const double       maxscale,
    const bool         printProgress,
    const Print      & print
) {

  PointContainer points = TdaToStl< PointContainer >(X, nSample, nDim,
      is_row_names);
  //lol copy paste
  //read_points(infilename, points);
  //read_points2(infilename, points);

  Distances distances(points); //PairDistances2
  Generator rips(distances); //Generator2
  typename Generator::Evaluator size(distances);
  Filtration filtration;
  EvaluatePushBack2< Filtration, typename Generator::Evaluator > functor(filtration, size);
  //auto functor = [&filtration](typename Generator::Simplex&& s) { filtration.push_back(s); };
  // Generate maxdimension skeleton of the Rips complex
  // rips.generate(skeleton, max_distance, [&filtration](Simplex&& s) { filtration.push_back(s); });
  //   
  rips.generate(maxdimension + 1, maxscale, functor);

  if (printProgress) {
    print("# Generated complex of size: %d \n", filtration.size());
  }

  // Sort the simplices with respect to comparison criteria
  // e.g. distance or function values
  // filtration.sort(ComparisonDataDimension< typename Filtration::Simplex >());
  filtration.sort(typename Generator::Comparison(distances));

  return filtration;
}


#endif __DIONYSUS2UTILS_H__
