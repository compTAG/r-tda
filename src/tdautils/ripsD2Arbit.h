//dionysus2
#include <dionysus/distances.h>
#include <geometry/Arbitdistance.h>
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

typedef         d::PairwiseDistances<std::vector<std::vector<double>>, ArbitDistance>        PairDistances2A;
typedef         PairDistances2A::DistanceType                             DistanceType2A;
typedef         PairDistances2A::IndexType                                VertexR2A;
typedef         d::Rips< PairDistances2A, d::Simplex< VertexR2A, double > >      Generator2A;
typedef         Generator2A::Simplex                                      SmplxR2A;
typedef         d::Filtration<SmplxR2A>                                       FltrR2A;
//typedef         StaticPersistence<>                                     PersistenceR;
//typedef         DynamicPersistenceChains<>                              PersistenceR;
//typedef         PersistenceDiagram<>                                    PDgmR;
typedef d::Z2Field K;
typedef d::PairRecorder<d::OrdinaryPersistence<K>> RipsPersistence2A;
