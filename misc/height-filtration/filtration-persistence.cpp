#include <iostream>
#include <unordered_map>

#include <dionysus/simplex.h>
#include <dionysus/filtration.h>
#include <dionysus/fields/z2.h>

#include "point.h"
#include "height_filtration.h"

typedef float Coordinate;
typedef grid::Point<Coordinate, 2> Point;
typedef Point Direction;
typedef dionysus::Z2Field K;
typedef dionysus::Simplex<Point> Simplex;
typedef dionysus::Filtration<Simplex> Filtration;

typedef HeightFiltrationFactory<Simplex, Direction, Filtration> HeightFiltFactory;

int main() {
    K k;

    Point p1({0,0});
    Point p2({2,2});
    Point p3({1,3});
    Point p4({4,1});
    Point p5({5,2});
    Point p6({6,1});

    Simplex t({p1, p2, p3});
    Simplex e1({p2, p4});
    Simplex e2({p4, p6});
    Simplex e3({p6, p5});
    Simplex e4({p5, p4});

    Direction d({-1,0});

    HeightFiltFactory hff(d);
    Filtration filtration = hff({ t, e1, e2, e3, e4 });

    for (auto& s : filtration) {
        std::cout << s << " " <<  filtration.index(s) << std::endl;
        for (auto sb : s.boundary(k)) {
            std::cout << "   " << sb.element() << " * "
                <<  sb.index() << " at "
                << filtration.index(sb.index()) << std::endl;
        }
    }
}
