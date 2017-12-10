#pragma once
#include "fft.hh"

namespace fpsm
{

namespace fft
{

    void init(int len)
    {
        buf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * len);
        fp = fftw_plan_dft_1d(len, buf, buf, FFTW_FORWARD,  FFTW_ESTIMATE);
        bp = fftw_plan_dft_1d(len, buf, buf, FFTW_BACKWARD, FFTW_ESTIMATE);
    }

    void cleanup()
    {
        fftw_destroy_plan(fp);
        fftw_destroy_plan(bp);
        fftw_free(buf);
    }

    template <class T>
    void transform_3d(int rank, grid const& g, std::vector<std::complex<T>>& a, direction const& d)
    {
    }
}

}

