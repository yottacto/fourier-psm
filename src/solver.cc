#pragma once
#include "solver.hh"
#include "fft/fft.hh"

namespace fpsm
{
    solver::solver(int ncd, int bnpcd)
        : g{ncd, bnpcd}, f{g.npc}, psi{g.npc}
    {
        MPI::Init();
        rank = MPI::COMM_WORLD.Get_rank();
        // f.reserve(g.npc);
        // psi.reserve(g.npc);
        fft::init(g.npcd);
    }

    solver::~solver()
    {
        MPI::Finalize();
        fft::cleanup();
    }

    void solver::init(func_type const& func)
    {
        for (auto i = 0; i < g.npc; i++)
            f[i] = func(g.get_point(rank, i));

        // TODO barrier
        fft::transform_3d(rank, g, f, fft::forward);
    }

    void solver::iterate(int num)
    {
        for (auto i = 0; i < num; i++) {
            format_rhs();
            normalize_lhs();
            // TODO barrier
            fft::transform_3d(rank, g, psi, fft::backward);
        }
    }

    void solver::format_rhs()
    {
        for (auto i = 0; i < g.npc; i++) {
            psi[i] *= std::norm(psi[i]);
            psi[i] = f[i] - psi[i];
        }
    }

    void solver::normalize_lhs()
    {
        // TODO reduce sum of (real, img) of id = 0
        std::complex<double> sum;
        for (auto i = 0; i < g.npc; i++) {
            psi[i] /= normalize_factor(g.get_index(rank, i));
            if (rank == 0 && i == 0) continue;
            sum += psi[i];
        }
    }

    double solver::normalize_factor(index const& p)
    {
    }

    void solver::print()
    {
    }

}

