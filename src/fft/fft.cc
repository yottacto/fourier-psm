#pragma once
#include <mpi.h>
#include "fft.hh"
#include "../utils/tools.hh"
#include "../utils/type.hh"

namespace fpsm
{

namespace fft
{
    namespace utils
    {
        int project(int rank, dimension const& d)
        {
            auto index = g.get_core_index(rank);
            if (d == dimension::x) return index.x;
            if (d == dimension::y) return index.y;
            // if (d == dimension::z)
            return index.z;
        }

        bool front_half(int phase, int rank, dimension const& d)
        {
            int id = project(rank, d);
            int base = g.bncd;
            return !(id & (1 << (base - phase + 1)));
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
            for (auto i = 1, len = g.ncd / 2; i < phase; i++, len /= 2)
                id %= len;
            return id;
        }

        std::complex<double> phase_factor(int phase, int rank, int local_id, dimension const& d)
        {
            auto project_id = project(rank, d);
            auto pha_id = phase_id(phase, project_id);

            // calc W_N^k
            auto N = g.npd / (1 << (phase - 1));
            auto k = pha_id * g.npcd + local_id;

            return {
                std::cos(-2. * pi * k / N),
                std::sin(-2. * pi * k / N)
            };
        }
    }

    void init(grid const& gr)
    {
        g = gr;
        auto len = g.npcd;
        buf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * len);
        fp = fftw_plan_dft_1d(len, buf, buf, FFTW_FORWARD,  FFTW_ESTIMATE);
        bp = fftw_plan_dft_1d(len, buf, buf, FFTW_BACKWARD, FFTW_ESTIMATE);
        send.resize(len);
        recv.resize(len);
    }

    void cleanup()
    {
        fftw_destroy_plan(fp);
        fftw_destroy_plan(bp);
        fftw_free(buf);
    }

    template <class T>
    void transform_3d(
        int rank,
        std::vector<std::complex<T>>& a,
        direction const& d)
    {
        if (d == backward)
            linear_transform_factor(rank, a, d);

        for (auto i = 0; i < 3; i++) {
            auto dime = static_cast<dimension>(i);
            transform_1d(rank, a, d, dime);
        }

        if (d == backward)
            linear_transform_factor(rank, a, d);
    }

    template <class T>
    void transform_1d(
        int rank,
        std::vector<std::complex<T>>& a,
        direction const& fft_d,
        dimension const& dim_d)
    {
        auto len = g.npcd;
        for (auto d1 = 0; d1 < len; d1++) {
            for (auto d2 = 0; d2 < len; d2++) {
                for (auto i = 0; i < len; i++) {
                    auto local_index = index_by_dimension(d1, d2, i, dim_d);
                    send[i] = a[g.get_local_id(local_index)];
                }

                // TODO barrier ?
                // init
                for (auto delta = g.ncd / 2, phase = 1; delta >= 1; delta /= 2, phase++) {
                    auto dest = utils::dest_id(phase, rank, dim_d);
                    MPI::COMM_WORLD.Sendrecv(
                        &send.front(), len, MPI::DOUBLE_COMPLEX, dest, 99,
                        &recv.front(), len, MPI::DOUBLE_COMPLEX, dest, 99
                    );

                    for (auto i = 0; i < len; i++) {
                        if (utils::front_half(i, rank, dim_d)) {
                            send[i] += recv[i];
                        } else {
                            send[i] -= recv[i];
                            send[i] *= utils::phase_factor(phase, rank, i, dim_d);
                        }
                    }
                }
            }
        }
    }

    template <class T>
    void linear_transform_factor(
        int rank,
        std::vector<std::complex<T>>& a,
        direction const& d)
    {
    }

}

}

