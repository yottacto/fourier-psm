#pragma once
#include <vector>
#include <complex>
#include <fftw3.h>
#include "../grid.hh"

namespace fpsm
{

namespace fft
{
    enum direction { forward, backward };

    fftw_complex* buf;
    fftw_plan fp;
    fftw_plan bp;

    void init(int len);
    void cleanup();

    template <class T>
    void transform_3d(
        int rank,
        grid const& g,
        std::vector<std::complex<T>>& a,
        direction const& d
    );
}

}

