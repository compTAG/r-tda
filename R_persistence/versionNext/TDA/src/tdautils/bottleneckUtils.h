// #include <utilities/types.h>
// #include <boost/archive/binary_iarchive.hpp>
// #include <boost/serialization/map.hpp>
#include <topology/persistence-diagram.h>

typedef PersistenceDiagram<>                    PDgmB;

template <typename RealVec>
void read_diagram(PDgmB& dgm, const RealVec & pointsInput, int pointsNumberInput)
{
	int idx;
	for (idx = 0; idx < pointsNumberInput ; ++idx)
	{
		dgm.push_back(PDgmB::Point(pointsInput[2*idx],pointsInput[2*idx+1]));
	}
}
