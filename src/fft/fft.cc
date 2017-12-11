#include <iostream>
// FIXME

#include "fft.hh"
#include "../utils/constant.hh"

namespace fpsm
{

namespace fft
{

    std::vector<std::complex<double>> send;
    std::vector<std::complex<double>> recv;

    fftw_complex* buf;
    fftw_plan plan[2];

    grid g;

    namespace utils
    {
        int core_project(int rank, dimension const& d)
        {
            auto index = g.get_core_index(rank);
            if (d == dimension::x) return index.x;
            if (d == dimension::y) return index.y;
            // if (d == dimension::z)
            return index.z;
        }

        int local_project(int local_id, dimension const& d)
        {
            auto index = g.get_local_index(local_id);
            if (d == dimension::x) return index.x;
            if (d == dimension::y) return index.y;
            // if (d == dimension::z)
            return index.z;
        }

        bool front_half(int phase, int rank, dimension const& d)
        {
            int id = core_project(rank, d);
            int base = g.bncd;
            return !(id & (1 << (base - phase)));
        }

        int next_delta_id(int rank, dimension const& d, int delta)
        {
            auto index = g.get_core_index(rank);
            if (d == dimension::x) index.x += delta;
            if (d == dimension::y) index.y += delta;
            if (d == dimension::z) index.z += delta;
            return g.get_core_rank(index);
        }

        int prev_delta_id(int rank, dimension const& d, int delta)
        {
            return next_delta_id(rank, d, -delta);
        }

        int dest_id(int phase, int rank, dimension const& d)
        {
            auto delta = g.ncd / (1 << phase);
            if (front_half(phase, rank, d))
                return next_delta_id(rank, d, delta);
            else
                return prev_delta_id(rank, d, delta);
        }

        int phase_id(int phase, int id)
        {
            for (auto i = 1, len = g.ncd / 2; i <= phase; i++, len /= 2)
                id %= len;
            return id;
        }

        std::complex<double> phase_factor(int phase, int rank, int local_pid, dimension const& d)
        {
            auto core_project_id = core_project(rank, d);
            auto pha_id = phase_id(phase, core_project_id);

            // calc W_N^k
            auto N = g.npd / (1 << (phase - 1));
            auto k = pha_id * g.npcd + local_pid; // local_project(local_id, d);

            return {
                std::cos(-2. * pi * k / N),
                std::sin(-2. * pi * k / N)
            };
        }

        int signed_parity(int x)
        {
            return x & 1 ? -1 : 1;
        }
    }

    void init(grid const& gr)
    {
        g = gr;
        auto len = g.npcd;
        buf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * len);
        plan[0] = fftw_plan_dft_1d(len, buf, buf, FFTW_FORWARD,  FFTW_ESTIMATE);
        plan[1] = fftw_plan_dft_1d(len, buf, buf, FFTW_BACKWARD, FFTW_ESTIMATE);
        send.resize(len);
        recv.resize(len);
    }

    void cleanup()
    {
        fftw_destroy_plan(plan[0]);
        fftw_destroy_plan(plan[1]);
        fftw_free(buf);
    }

}

}

