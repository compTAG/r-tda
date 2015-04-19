// #include <utilities/types.h>
// #include <boost/archive/binary_iarchive.hpp>
// #include <boost/serialization/map.hpp>
#include <topology/persistence-diagram.h>

typedef PersistenceDiagram<>                    PDgmB;

template <typename RealMatrix>
inline PDgmB read_diagram(const RealMatrix& Diag)
{
	const unsigned diagNum = Diag.nrow();
	PDgmB dgm;
	for (unsigned diagIdx = 0; diagIdx < diagNum; ++diagIdx)
	{
		dgm.push_back(PDgmB::Point(Diag[diagIdx + 0 * diagNum], Diag[diagIdx + 1 * diagNum]));
	}
	return dgm;
}
