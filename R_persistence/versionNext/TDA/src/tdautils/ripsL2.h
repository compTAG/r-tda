#include <topology/rips.h>
#include <topology/filtration.h>
#include <topology/static-persistence.h>
#include <topology/dynamic-persistence.h>
#include <topology/persistence-diagram.h>

#include <geometry/l2distance.h>
#include <geometry/distances.h>
#include <utilities/containers.h>           // for BackInsertFunctor
#include <utilities/timer.h>

#include <vector>

typedef         PairwiseDistances<PointContainer, L2Distance>           PairDistances;
typedef         PairDistances::DistanceType                             DistanceType;
typedef         PairDistances::IndexType                                VertexR;
typedef         Rips<PairDistances>                                     Generator;
typedef         Generator::Simplex                                      SmplxR;
typedef         Filtration<SmplxR>                                       FltrR;
typedef         StaticPersistence<>                                     PersistenceR;
//typedef         DynamicPersistenceChains<>                              PersistenceR;
typedef         PersistenceDiagram<>                                    PDgmR;


