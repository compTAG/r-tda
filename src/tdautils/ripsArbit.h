#include <topology/rips.h>
//#include <topology/filtration.h>
#include <topology/static-persistence.h>
#include <topology/dynamic-persistence.h>
#include <topology/persistence-diagram.h>

//#include <geometry/l2distance.h>
#include <geometry/Arbitdistance.h>
//#include <geometry/distances.h>
#include <utilities/containers.h>           // for BackInsertFunctor
#include <utilities/timer.h>

//dionysus2
#include <dionysus/filtration.h>
#include <dionysus/distances.h>

#include <vector>

namespace d = dionysus;

typedef         d::PairwiseDistances<std::vector<std::vector<double>>, ArbitDistance>        PairDistancesA;
typedef         PairDistancesA::DistanceType                             DistanceTypeA;
typedef         PairDistancesA::IndexType                                VertexRA;
typedef         Rips< PairDistancesA, Simplex< VertexRA, double > >      GeneratorA;
typedef         GeneratorA::Simplex                                      SmplxRA;
typedef         Filtration<SmplxRA>                                       FltrRA;
typedef         StaticPersistence<>                                     PersistenceR;
//typedef         DynamicPersistenceChains<>                              PersistenceR;
typedef         PersistenceDiagram<>                                    PDgmR;
