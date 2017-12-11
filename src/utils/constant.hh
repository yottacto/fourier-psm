#pragma once
#include <boost/math/constants/constants.hpp>
#include <cmath>
#include "type.hh"

namespace fpsm
{
    auto constexpr pi = boost::math::constants::pi<double>();

    inline std::complex<double> default_f(point const& p)
    {

        return {
            -0.5 + 1.5 * cos(p.z) - 3. * cos(p.y) + 2. * cos(p.x),
            1.5 * sin(p.z) - 3. * sin(p.y) + 2. * sin(p.x),
        };

        return {
            std::cos(p.x) * std::cos(p.y) * std::cos(p.z),
            std::sin(p.x) * std::sin(p.y) * std::sin(p.z)
        };
    };
}

