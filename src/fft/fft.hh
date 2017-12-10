#pragma once
#include <vector>
#include <complex>
#include <cmath>
#include <fftw3.h>
#include "../grid.hh"
#include "../utils/type.hh"

namespace fpsm
{

namespace fft
{
    enum direction { forward, backward };

    std::vector<std::complex<double>> send;
    std::vector<std::complex<double>> recv;

    fftw_complex* buf;
    fftw_plan plan[2];

    grid g;

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
    }

    template <class T>
    void linear_transform_factor(
        int rank,
        std::vector<std::complex<T>>& a
    );

    template <class T>
    void transform_1d(
        int rank,
        std::vector<std::complex<T>>& a,
        direction const& fft_d,
        dimension const& dim_d
    );

    template <class T>
    void transform_3d(
        int rank,
        std::vector<std::complex<T>>& a,
        direction const& d
    );
}

}

