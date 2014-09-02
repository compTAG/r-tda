#include <utilities/types.h>
#include <string>
#include <sstream>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/map.hpp>


typedef PersistenceDiagram<>                    PDgmB;

void read_diagram(PDgmB& dgm, const std::string& filename)
{
    std::ifstream in(filename.c_str());
    std::string line;
    std::getline(in, line);
    while (in)
    {
        std::istringstream sin(line);
        double x,y;
        sin >> x >> y;
        dgm.push_back(PDgmB::Point(x,y));
        std::getline(in, line);
    }
}
