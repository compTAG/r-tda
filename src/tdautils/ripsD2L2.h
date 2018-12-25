//dionysus2
#include <dionysus/distances.h>
#include <dionysus/fields/z2.h>
#include <dionysus/filtration.h>
#include <dionysus/rips.h>
#include <dionysus/ordinary-persistence.h>
#include <dionysus/standard-reduction.h>
#include <dionysus/diagram.h>
#include <dionysus/pair-recorder.h>

#include <utilities/containers.h>           // for BackInsertFunctor
#include <utilities/timer.h>


#include <vector>

namespace d = dionysus;


//L2 Struct is inside distances.h

typedef     std::vector<double>                                     Point;
typedef     std::vector<Point>                                      PointContainer;

typedef         d::PairwiseDistances<PointContainer, d::L2Distance<Point>>   PairDistances2;
typedef         PairDistances2::DistanceType                             DistanceType2;
typedef         PairDistances2::IndexType                                VertexR2;

typedef         d::Rips< PairDistances2, d::Simplex< VertexR, double > >       Generator2;
typedef         Generator2::Simplex                                      SmplxR2;
typedef         d::Filtration<SmplxR2>                                      FltrR2;

typedef d::Z2Field K;
typedef d::PairRecorder<d::OrdinaryPersistence<K>> RipsPersistence2;
