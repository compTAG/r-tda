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



/*
     * I'm going to assume that we simply pass in the persDgm, persLoc, and persCycle as they are above.
     * I should ask Dave tomorrow what persDgm, persLoc, and persCycle are.
     * Persistence Locations: for each diagram, there is a collection of points, and the points are "steps" 
     * Perstistence Cycle: Birth and Death?
     *
     * Make example to test it
     *
     * ask how to get the vignette
     * Figure out what the Filtration that gets passed in is.
     * Want to pass in Ordinary-Persistence with a Z2Field first go
     * 
     * ask dave what initLocations and initDiagrams do
     *
     * ask dave what format the cmplx is passed in as
     * Figure out what format Ordinary Persistence saves the persistence as
     *
     * Figure out what format TDA filtration is in typecast utils
     * What is the format for TDA filtration?
     *
     * Figure out what <Persistence Parameter> is in FiltrationDiagDionysus
     *      Parameter probably comes from gridUtils.h defining Persistence as static-persistence
     *      Probably want to rename persistence2
     * Figure out how to convert Filtration.h in D1 to D2 so the rest of the methods work.
     * 
     */
    
//Helper function for filling in persDgm
//template< typename Diagrams, typename iterator, typename Evaluator, typename SimplexMap >
//inline void initDiagrams;
//Helper function for filling in persLoc and persCycle
// inline void initLocations;

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
    
    //Assume that Persistence that is passed in is Persistence2
    //Calculate Persistence 
    
    d::Z2Field k;
    //Persistence persistence(k);
    //StandardReduction2 reduce(persistence);
    d::RowReduction<d::Z2Field> reduce(k);
    // We know that the function breaks when this line is called.
    reduce(filtration);
    
    typedef decltype(reduce.persistence().pair(0)) Index;
    typedef float Value;
    //persistence is reduced.
    Index _ = 0;
    // move Persistence into persDgm
    //auto dgms = d::init_diagrams(reduce->persistence(), filtration, [&](const Smplx2& s) -> float  { return filtration.index(s); }, [](typename Persistence::Index i) { return i; });
    
    diagramDS::DiagramDS<float, Index> dgms(
            reduce.persistence(), 
            filtration, 
            [&](const Smplx2& s) -> float  { return filtration.index(s);}, 
            [](typename Persistence::Index i) { return i; }
        );
    //emulate initDiagrams function from dionysusUtils
    //will put into a function later
    persDgm.resize(dgms.getDiagrams().size());
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
    persDgm.resize(maxdimension + 1);
    
}

#endif __DIONYSUS2UTILS_H__
