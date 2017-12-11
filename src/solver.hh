#pragma once
#include <vector>
#include <complex>
#include <mpi.h>
#include "grid.hh"
#include "utils/type.hh"

namespace fpsm
{
    struct solver
    {
        solver(int ncd, int bnpcd);
        ~solver();

        void init(func_type const& func);
        void iterate(int num);
        void print() const;
        void print(std::vector<std::complex<double>> const& a) const;

        void normalize_lhs();
        void format_rhs();
        double normalize_factor(index p) const;

        int rank;
        grid g;
        std::vector<std::complex<double>> f;
        std::vector<std::complex<double>> psi;
    };
}

