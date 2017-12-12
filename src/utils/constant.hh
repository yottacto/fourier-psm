#pragma once
#include <boost/math/constants/constants.hpp>
#include <cmath>
#include "type.hh"

namespace fpsm
{
    auto constexpr pi = boost::math::constants::pi<double>();

    inline std::complex<double> debug_f(point const& p)
    {
        return {
            p.x + p.y + p.z,
            p.x * p.y * p.z
        };
    }

    inline std::complex<double> test_f(point const& p)
    {
        return {
            -0.5 + 1.5 * std::cos(p.z) - 3. * std::cos(p.y) + 2. * std::cos(p.x),
            1.5        * std::sin(p.z) - 3. * std::sin(p.y) + 2. * std::sin(p.x),
        };
    }

    inline std::complex<double> default_f(point const& p)
    {
        return {
            std::sin(p.x) * std::sin(p.y) * std::sin(p.z),
            std::cos(p.x) * std::cos(p.y) * std::cos(p.z)
        };
    };
}

