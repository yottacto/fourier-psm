#pragma once
#include <functional>
#include <type_traits>
#include <complex>

namespace fpsm
{
    struct point
    {
        double x;
        double y;
        double z;
    };

    struct index
    {
        int x;
        int y;
        int z;
    };

    using func_type = std::function<std::complex<double>(point const&)>;
    using func_result_type = std::result_of_t<func_type(point)>;
}

