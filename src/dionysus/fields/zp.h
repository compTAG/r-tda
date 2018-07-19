#ifndef DIONYSUS_ZP_H
#define DIONYSUS_ZP_H

#include <vector>

namespace dionysus
{

template<typename Element_ = short>
class ZpField
{
    public:
        typedef         Element_                            Element;
        
                        ZpField(Element p);
                        ZpField(const ZpField& other)       = default;
                        ZpField(ZpField&& other)            = default;

        Element         id()  const                         { return 1; }
        Element         zero()  const                       { return 0; }
        // adding a prime number so that when modded the result is positive
        Element         init(int a) const                   { return (a % p_ + p_) % p_; }

        Element         neg(Element a) const                { return p_ - a; }
        Element         add(Element a, Element b) const     { return (a+b) % p_; }

        Element         inv(Element a) const                { 
            while (a < 0) a += p_;
            return inverses_[a]; 
        }
        
        Element         mul(Element a, Element b) const     { return (a*b) % p_; }
        Element         div(Element a, Element b) const     { return mul(a, inv(b)); }

        bool            is_zero(Element a) const            { return (a % p_) == 0; }

        Element         prime() const                       { return p_; }

    private:
        Element                 p_;
        std::vector<Element>    inverses_;
};

template<class E>
ZpField<E>::ZpField(Element p): 
    p_(p), // constructor for setting p_
    inverses_(p_) // constructor that sets length of vector
    {
    for (Element i = 1; i < p_; ++i)
        for (Element j = 1; j < p_; ++j)
            if (mul(i,j) == 1)
            {
                inverses_[i] = j;
                break;
            }
    }

}

#endif
