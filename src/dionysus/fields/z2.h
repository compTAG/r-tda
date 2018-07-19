#ifndef DIONYSUS_Z2_H
#define DIONYSUS_Z2_H

namespace dionysus
{

class Z2Field
{
    public:
        typedef         short                               Element; 

                        Z2Field()                           {} // this is a constructor
        // this is a function that returns a short 1
        static Element  id()                                { return 1; } 
        // this is a function that returns 0
        static Element  zero()                              { return 0; }         
        // this init function returns 
        static Element  init(int a)                         { return (a % 2 + 2) % 2; } 
        // turn a 0 to a 1 or vice versa
        Element         neg(Element a) const                { return 2 - a; }
        // add elements in binary
        Element         add(Element a, Element b) const     { return (a+b) % 2; }
        Element         inv(Element a) const                { return a; }
        Element         mul(Element a, Element b) const     { return a*b; }
        // This is strange
        Element         div(Element a, Element b) const     { return a; }

        bool            is_zero(Element a) const            { return a == 0; }
};

}

#endif

