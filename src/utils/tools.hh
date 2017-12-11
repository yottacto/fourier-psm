#pragma once
#include "type.hh"

namespace fpsm
{
    inline index index_by_dimension(int base1, int base2, int d, dimension const& dime)
    {
        if (dime == dimension::x) return {d, base1, base2};
        if (dime == dimension::y) return {base1, d, base2};
        // if (dime == dimension::z)
        return {base1, base2, d};
    }
}

