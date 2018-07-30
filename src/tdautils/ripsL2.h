#include <topology/rips.h>
//#include <topology/filtration.h>
#include <topology/static-persistence.h>
#include <topology/dynamic-persistence.h>
#include <topology/persistence-diagram.h>

//#include <geometry/l2distance.h>
//#include <geometry/distances.h>
#include <utilities/containers.h>           // for BackInsertFunctor
#include <utilities/timer.h>
 
//dionysus2
//#include <dionysus/rips.h>
#include <dionysus/filtration.h>
//#include <dionysus/ordinary-persistence.h>
#include <dionysus/distances.h>
//#include <dionysus/diagram.h>

#include <vector>

namespace d = dionysus;


//L2 Struct is inside distances.h
// not sure if I want to typedef here

// typedef     std::vector<double>                                     Point;
// typedef     std::vector<Point>                                      PointContainer;

// Swapped out
typedef         d::PairwiseDistances<std::vector<std::vector<double>>, d::L2Distance<std::vector<double>>>           PairDistances;

//Next two lines are fine
typedef         PairDistances::DistanceType                             DistanceType;
typedef         PairDistances::IndexType                                VertexR;


typedef         Rips< PairDistances, Simplex< VertexR, double > >       Generator;
typedef         Generator::Simplex                                      SmplxR;
typedef         Filtration<SmplxR>                                       FltrR;

// Comment this test
typedef         StaticPersistence<>                                     PersistenceR;

// relabel
//typedef         OrdinaryPersistence<>                                   PersistenceR;

//typedef         DynamicPersistenceChains<>                              PersistenceR;
typedef         PersistenceDiagram<>                                    PDgmR;
 
