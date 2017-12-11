#pragma once
#include <iomanip>
#include <fstream>
#include "solver.hh"
#include "fft/fft.hh"

namespace fpsm
{
    solver::solver(int ncd, int bnpcd)
        : g{ncd, bnpcd}, f(g.npc), psi(g.npc)
    {
        MPI::Init();
        rank = MPI::COMM_WORLD.Get_rank();
        // f.reserve(g.npc);
        // psi.reserve(g.npc);
        fft::init(g);
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
        // fft::transform_3d(rank, f, fft::forward);

        print(f);
    }

    void solver::iterate(int num)
    {
        for (auto i = 0; i < num; i++) {
            format_rhs();
            normalize_lhs();
            // TODO barrier
            fft::transform_3d(rank, psi, fft::backward);
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

    double solver::normalize_factor(index p) const
    {
        if (p.x >= g.npd/2) p.x -= g.npd;
        if (p.y >= g.npd/2) p.y -= g.npd;
        if (p.z >= g.npd/2) p.z -= g.npd;
        return p.x * p.x + p.y * p.y + p.z * p.z;
    }

    void solver::print(std::vector<std::complex<double>> const& a) const
    {
        std::ofstream fout{"output/out.dat"};
        std::vector<std::complex<double>> out;
        if (rank == 0) out.resize(g.np);
        MPI::COMM_WORLD.Gather(
            &a.front(),   g.npc, MPI::DOUBLE_COMPLEX,
            &out.front(), g.npc, MPI::DOUBLE_COMPLEX,
            0
        );
        if (rank == 0) {
            auto tmp = out;
            for (auto i = 0; i < g.nc; i++) {
                for (auto lid = 0; lid < g.npc; lid++) {
                    auto id = g.get_id(i, lid);
                    out[id] = tmp[i * g.npc + lid];
                }
            }

            fout << "n = " << g.np << "\n";
            for (auto i = 0; i < g.np; i++) {
                auto p = g.get_point(i);
                fout << p.x << " " << p.y << " " << p.z << " "
                    << out[i].real() << " " << out[i].imag() << "\n";
            }
        }

    }

    void solver::print() const
    {
        print(psi);
    }

}

