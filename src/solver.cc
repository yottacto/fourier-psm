#include <iomanip>
#include <fstream>
#include "solver.hh"
#include "fft/fft.hh"

namespace fpsm
{
    solver::solver(int ncd, int bnpcd)
        : g{ncd, bnpcd}, f(g.npc), psi(g.npc)
    {
        rank = MPI::COMM_WORLD.Get_rank();
        // f.reserve(g.npc);
        // psi.reserve(g.npc);
        fft::init(g);
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

            // // TODO remove
            // auto p = g.get_point(rank, i);
            // if (rank == 0 && i == 0) {
            //     std::cout << p.x << ", " << p.y << ", " << p.z << "\n";
            //     std::cout << f[i] << "\n";
            // }

        }

        // TODO remove this
        // if (rank == 0) {
        //     for (int i = 0; i < g.nc; i++)
        //         std::cerr << i << " -> " <<
        //             fft::utils::dest_id(1, i, dimension::y) << "\n";
        // }
        // return;

        fft::transform_3d(rank, f, fft::forward);

        // print(f);
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
            psi[i] = f[i] - psi[i];
        }
    }

    void solver::normalize_lhs()
    {
        // TODO reduce sum of (real, img) of id = 0
        std::complex<double> sum;
        for (auto i = 0; i < g.npc; i++) {
            if (rank == 0 && i == 0) continue;
            auto nf = normalize_factor(g.get_index(rank, i));
            psi[i] /= nf;
            sum += psi[i];
        }
        std::complex<double> factor;
        // FIXME barrier?
        MPI::COMM_WORLD.Barrier();
        MPI::COMM_WORLD.Reduce(
            &sum, &factor, 1, MPI::DOUBLE_COMPLEX, MPI::SUM, 0
        );
        if (rank == 0) {
            psi[0] = -factor;
        }
    }

    double solver::normalize_factor(index p) const
    {
        if (p.x >= g.npd/2) p.x -= g.npd;
        if (p.y >= g.npd/2) p.y -= g.npd;
        if (p.z >= g.npd/2) p.z -= g.npd;
        return p.x * p.x + p.y * p.y + p.z * p.z;
    }

    // TODO remove
    template <class T>
    T abs(T a) { return a < 0 ? -a : a; }

    template <class T>
    T knight(T a)
    {
        auto const eps = 1e-8;
        if (abs(a.real()) < eps) a.real(0);
        if (abs(a.imag()) < eps) a.imag(0);
        return a;
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

            // // TODO remove
            // std::cerr << knight(out[g.get_id({0, 0, 0})]) << "\n";
            // std::cerr << knight(out[g.get_id({0, 0, 1})]) << "\n";
            // std::cerr << knight(out[g.get_id({0, 1, 0})]) << "\n";
            // std::cerr << knight(out[g.get_id({1, 0, 0})]) << "\n";
            // std::cerr << knight(out[g.get_id({20, 23, 22})]) << "\n";
            // std::cerr << knight(out[g.get_id({17, 23, 11})]) << "\n";

            for (int i = 0; i < g.np; i++)
                std::cerr << knight(out[i]) << "\n";

            for (int i = 0; i < 20; i++)
                std::cerr << knight(out[i]) << "\n";

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
        print(f);
    }

}

