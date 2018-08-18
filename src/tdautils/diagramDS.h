#ifndef DIAGRAM_DS_H
#define DIAGRAM_DS_H
#endif

#include <vector>
#include <tuple>

namespace d = dionysus;

namespace diagramDS
{

template<class Value_, class Data_>
class DiagramDS
{
    public:
        using Value = Value_;
        using Data = Data_; 
        using Diagrams = std::vector<d::Diagram<Value,Data>>;
        
        template<class ReducedMatrix, class Filtration, class GetValue, class GetData>
        DiagramDS(const ReducedMatrix& m, const Filtration& f, const GetValue get_value, const GetData get_data)
        {
            //auto get_value = [&](const Simplex& s) -> float  { return filtration.index(s); };
            //auto get_data = [](Persistence::Index i) { return i; };
            for (typename ReducedMatrix::Index i = 0; i < m.size(); ++i)
            {
                if (m.skip(i))
                    continue;

                auto& s = f[i];
                auto d = s.dimension();
                while (d + 1 > diagrams.size())
                    diagrams.emplace_back();

                auto pair = m.pair(i);
                if (pair == m.unpaired())
                {
                    auto birth = get_value(s);
                    using Value = decltype(birth);
                    Value death = std::numeric_limits<Value>::infinity();
                    diagrams[d].emplace_back(birth, death, get_data(i));
                } else if (pair > i)       // positive
                {
                    auto birth = get_value(s);
                    auto death = get_value(f[pair]);

                    if (birth != death)         // skip diagonal
                        diagrams[d].emplace_back(birth, death, get_data(i));
                } // else negative: do nothing
            }
        }
        
        Diagrams getDiagrams() {return diagrams;}
    private:
        Diagrams diagrams;
};

}


