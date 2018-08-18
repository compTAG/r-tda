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
            persDgm[_].push_back(pt_);
        }
        _++;   
    }
    // Capping at maxdimension
    if (persDgm.size() > maxdimension) {
        persDgm.resize(maxdimension + 1);
    }
}

#endif __DIONYSUS2UTILS_H__
