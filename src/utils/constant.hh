#pragma once
#include <boost/math/constants/constants.hpp>
#include <cmath>
#include "type.hh"

namespace fpsm
{
    auto constexpr pi = boost::math::constants::pi<double>();

    auto constexpr default_f = [](point const& p) {
        return std::complex<double>{
            std::cos(p.x) * std::cos(p.y) * std::cos(p.z),
            std::sin(p.x) * std::sin(p.y) * std::sin(p.z)
        };
    };
}

