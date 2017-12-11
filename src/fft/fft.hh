#pragma once
#include <vector>
#include <complex>
#include <cmath>
#include <fftw3.h>
#include <mpi.h>
#include "../grid.hh"
#include "../utils/type.hh"
#include "../utils/tools.hh"

namespace fpsm
{

namespace fft
{
    enum direction { forward, backward };

    extern std::vector<std::complex<double>> send;
    extern std::vector<std::complex<double>> recv;

    extern fftw_complex* buf;
    extern fftw_plan plan[2];

    extern grid g;

    void init(grid const& g);
    void cleanup();

    namespace utils
    {
        int project(int rank, dimension const& d);
        int dest_id(int phase, int rank, dimension const& d);
        bool front_half(int phase, int rank, dimension const& d);
        std::complex<double> phase_factor(int phase, int rank, int local_id, dimension const& d);
        int next_delta_id(int rank, dimension const& d, int delta);
        int prev_delta_id(int rank, dimension const& d, int delta);
        int phase_id(int phase, int id);
        int signed_parity(int x);
    }

    template <class T>
    void linear_transform_factor(int rank, std::vector<std::complex<T>>& a)
    {
        for (auto i = 0; i < g.npc; i++) {
            auto index = g.get_index(rank, i);
            a[i] *= utils::signed_parity(index.x + index.y + index.z);
        }
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

                // init
                auto phase = 1;
                for (auto delta = g.ncd / 2; delta >= 1; delta /= 2, phase++) {

                    auto dest = utils::dest_id(phase, rank, dim_d);

                    MPI::COMM_WORLD.Barrier();
                    MPI::COMM_WORLD.Sendrecv(
                        &send.front(), len, MPI::DOUBLE_COMPLEX, dest, 99,
                        &recv.front(), len, MPI::DOUBLE_COMPLEX, dest, 99
                    );

                    if (utils::front_half(phase, rank, dim_d)) {
                        for (auto i = 0; i < len; i++)
                            send[i] += recv[i];
                    } else {
                        for (auto i = 0; i < len; i++) {
                            send[i] -= recv[i];
                            send[i] *= utils::phase_factor(phase, rank, i, dim_d);
                        }
                    }
                }

                // main fftw
                for (auto i = 0; i < len; i++) {
                    buf[i][0] = recv[i].real();
                    buf[i][1] = recv[i].imag();
                }

                fftw_execute(plan[fft_d]);

                // reindex, reverse phase
                for (auto i = 0; i < len; i++)
                    send[i] = {buf[i][0], buf[i][1]};

                phase--;
                for (auto delta = 1; delta <= g.ncd/2; delta *= 2, phase--) {
                    auto dest = utils::dest_id(phase, rank, dim_d);

                    // TODO remove this debug
                    // if (rank >= 60) {
                    //     std::cerr << "phase = " << phase <<
                    //         " rank = " << rank <<
                    //         " dest = " << dest << std::endl;
                    //     // std::cerr << "[rank  == " << rank << "] [phase == " << phase << "] [dest == " << dest << "]" << std::endl;
                    // }

                    MPI::COMM_WORLD.Barrier();
                    MPI::COMM_WORLD.Sendrecv(
                        &send.front(), len, MPI::DOUBLE_COMPLEX, dest, 99,
                        &recv.front(), len, MPI::DOUBLE_COMPLEX, dest, 99
                    );
                    if (utils::front_half(phase, rank, dim_d)) {
                        for (auto i = 0; i < len/2; i++)
                            recv[i + len/2] = send[i];
                        for (auto i = 0; i < len/2; i++) {
                            send[i * 2]     = recv[i + len/2];
                            send[i * 2 + 1] = recv[i];
                        }
                    } else {
                        for (auto i = 0; i < len/2; i++)
                            recv[i] = send[i + len/2];
                        for (auto i = 0; i < len/2; i++) {
                            send[i * 2]     = recv[i + len/2];
                            send[i * 2 + 1] = recv[i];
                        }
                    }
                }

                // place back to a
                for (auto i = 0; i < len; i++) {
                    auto local_index = index_by_dimension(d1, d2, i, dim_d);
                    // FIXME (T)
                    auto scale = 1.;
                    if (fft_d == forward) scale = g.npd;
                    a[g.get_local_id(local_index)] = send[i] / scale;
                }
            }
        }
    }

    template <class T>
    void transform_3d(
        int rank,
        std::vector<std::complex<T>>& a,
        direction const& d)
    {
        if (d == backward)
            linear_transform_factor(rank, a);

        for (auto i = 0; i < 3; i++) {
            auto dime = dimension(i);
            transform_1d(rank, a, d, dime);
        }

        if (d == forward)
            linear_transform_factor(rank, a);
    }
}

}

