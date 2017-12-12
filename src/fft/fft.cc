#include "fft.hh"
#include "../utils/constant.hh"

namespace fpsm
{

namespace fft
{

    std::vector<std::complex<double>> send;
    std::vector<std::complex<double>> recv;

    std::vector<std::complex<double>> reorder;
    std::vector<std::complex<double>> tmp_buf;
    std::vector<int> index;

    fftw_complex* buf;
    fftw_plan plan[2];

    grid g;

    namespace utils
    {
        int head_of_core_in_dimension(int rank, dimension const& d)
        {
            auto index = g.get_core_index(rank);
            if (d == dimension::x) index.x = 0;
            if (d == dimension::y) index.y = 0;
            if (d == dimension::z) index.z = 0;
            return g.get_core_rank(index);
        }

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

        std::complex<double> phase_factor(int phase, int rank, int local_pid, dimension const& d, direction const& fft_d)
        {
            auto core_project_id = core_project(rank, d);
            auto pha_id = phase_id(phase, core_project_id);

            // calc W_N^k
            auto N = double(g.npd) / (1 << (phase - 1));
            auto k = double(pha_id * g.npcd + local_pid);

            auto scale = double{1};
            if (fft_d == direction::forward)
                scale = -1.;

            return {
                std::cos(scale * 2. * pi * k / N),
                std::sin(scale * 2. * pi * k / N)
            };
        }

        int signed_parity(int x)
        {
            return x & 1 ? -1 : 1;
        }
    }

    namespace impl
    {
        int left_circle_shift(int x, int len)
        {
            auto low = 0;
            if (x & (1 << (len - 1))) low = 1;
            x ^= low << (len - 1);
            x <<= 1;
            x |= low;
            return x;
        }

        void calc_index()
        {
            for (auto i = 0; i < g.npd; i++)
                index[i] = i;
            auto len = g.bncd + g.bnpcd;
            for (auto phase = 0, d = g.ncd / 2; d >= 1; d /= 2, phase++) {
                for (auto i = 0; i < g.npd; i++) {
                    auto base = index[i] & ((1 << phase) - 1);
                    index[i] >>= phase;
                    index[i] = left_circle_shift(index[i], len - phase);
                    index[i] <<= phase;
                    index[i] |= base;
                }
            }
        }

        bool can_be_head_core(int rank)
        {
            auto index = g.get_core_index(rank);
            return index.x == 0 || index.y == 0 || index.z == 0;
        }
    }

    void init(int rank, grid const& gr)
    {
        g = gr;
        auto len = g.npcd;
        buf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * len);
        plan[0] = fftw_plan_dft_1d(len, buf, buf, FFTW_FORWARD,  FFTW_ESTIMATE);
        plan[1] = fftw_plan_dft_1d(len, buf, buf, FFTW_BACKWARD, FFTW_ESTIMATE);
        send.resize(len);
        recv.resize(len);

        if (impl::can_be_head_core(rank)) {
            tmp_buf.resize(g.npd);
            reorder.resize(g.npd);
            index.resize(g.npd);
            impl::calc_index();
        } else {
            tmp_buf.resize(1);
            reorder.resize(1);
        }
    }

    void cleanup()
    {
        fftw_destroy_plan(plan[0]);
        fftw_destroy_plan(plan[1]);
        fftw_free(buf);
    }
}

}

