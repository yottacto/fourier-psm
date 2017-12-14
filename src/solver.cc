#include <iomanip>
#include <fstream>
#include "solver.hh"
#include "fft/fft.hh"

namespace fpsm
{
    solver::solver(int ncd, int bnpcd)
        : g{ncd, bnpcd}, f(g.npc), psi(g.npc), psi_transformed(g.npc)
    {
        rank = MPI::COMM_WORLD.Get_rank();
        fft::init(rank, g);
    }

    solver::~solver()
    {
        fft::cleanup();
    }

    void solver::init(func_type const& func)
    {
        for (auto i = 0; i < g.npc; i++) {
            f[i] = func(g.get_point(rank, i));
            psi[i] = {0, 0};
            psi_transformed[i] = {0, 0};
        }

        fft::transform_3d(rank, f, fft::forward);
        fft::transform_3d(rank, psi_transformed, fft::forward);
    }

    void solver::iterate(int num)
    {
        for (auto i = 0; i < num; i++) {
            format_rhs();
            normalize_lhs();
            fft::transform_3d(rank, psi, fft::backward);
        }
    }

    void solver::format_rhs()
    {
        for (auto i = 0; i < g.npc; i++) {
            psi[i] *= std::norm(psi[i]);
        }

        fft::transform_3d(rank, psi, fft::forward);

        for (auto i = 0; i < g.npc; i++) {
            psi[i] = f[i] - psi[i] + psi_transformed[i];
        }
    }

    void solver::normalize_lhs()
    {
        for (auto i = 0; i < g.npc; i++) {
            // if (rank == 0 && i == 0) continue;
            auto nf = normalize_factor(g.get_index(rank, i)) + 1.;
            psi[i] /= nf;
        }
        psi_transformed = psi;
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

        if (rank == 0)
            out.resize(g.np);
        // in case of reference binding to null pointer
        else
            out.resize(1);

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

            fout << "Title=\"solution\"\n";
            fout << "Variables=\"x\", \"y\", \"z\", \"real\", \"imag\"\n";
            fout << "ZONE I=" << g.npd << " J=" << g.npd << " K=" << g.npd << " f=point\n";
            for (auto i = 0; i < g.np; i++) {
                auto p = g.get_point(i);
                fout << p.x << " " << p.y << " " << p.z << " "
                    << out[i].real() << " " << out[i].imag() << "\n";
            }
        }

    }

    void solver::print(func_type const& func) const
    {
        std::ofstream fout{"output/out_func.dat"};

        if (rank == 0) {
            fout << "Title=\"solution\"\n";
            fout << "Variables=\"x\", \"y\", \"z\", \"real\", \"imag\"\n";
            fout << "ZONE I=" << g.npd << " J=" << g.npd << " K=" << g.npd << " f=point\n";

            for (auto i = 0; i < g.np; i++) {
                auto p = g.get_point(i);
                auto value = func(p);
                fout << p.x << " " << p.y << " " << p.z << " "
                    << value.real() << " " << value.imag() << "\n";
            }
        }
    }

    void solver::print() const
    {
        print(psi);
    }

}

