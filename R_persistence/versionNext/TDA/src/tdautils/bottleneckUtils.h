// #include <utilities/types.h>
// #include <boost/archive/binary_iarchive.hpp>
// #include <boost/serialization/map.hpp>
#include <topology/persistence-diagram.h>

typedef PersistenceDiagram<>                    PDgmB;

template <typename RealMatrix>
void read_diagram(PDgmB& dgm, const RealMatrix& Diag)
{
	int idx;
	for (idx = 0; idx < Diag.nrow(); ++idx)
	{
		dgm.push_back(PDgmB::Point(Diag(idx, 0), Diag(idx, 1)));
	}
}
