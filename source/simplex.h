#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <boost/dynamic_bitset.hpp>

using namespace boost;

class implicitS
{
    int top_dim;
    int top_index;
    dynamic_bitset<> simpl;

public:
    Simplex(int,int,int);
};

#endif // SIMPLEX_H
