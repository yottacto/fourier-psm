#include <iostream>
#include <mpi.h>
#include "solver.hh"
#include "utils/constant.hh"
#include "fft/fft.hh"
#include "utils/type.hh"
#include "config.hh"

int main()
{
    MPI::Init();
    auto rank = MPI::COMM_WORLD.Get_rank();

    auto toml = cpptoml::parse_file("config.toml");
    icesp::configuration::config config{toml};
    // base number of cores per direction
    int bncd{config.cases.begin()->second.bncd};
    // base number of points per core per direction
    int bnpcd{config.cases.begin()->second.bnpcd};
    int iteration{config.cases.begin()->second.iteration};
    fpsm::solver s(bncd, bnpcd);
    if (rank == 0) {
        std::cout << "n = " << s.g.npd << "\n";
        std::cout << "tot = " << s.g.np << "\n";
    }

    MPI::COMM_WORLD.Barrier();
    auto start = MPI::Wtime();

    // s.init(fpsm::debug_f);
    s.init(fpsm::default_f);
    s.iterate(iteration);
    s.print();
    // s.print(fpsm::default_f);

    MPI::COMM_WORLD.Barrier();
    auto end = MPI::Wtime();
    auto elapsed = end - start;
    auto tmp = elapsed;
    MPI::COMM_WORLD.Reduce(&tmp, &elapsed, 1, MPI::DOUBLE, MPI::MAX, 0);
    MPI::Finalize();

    if (rank == 0)
        std::cout << "elapsed time: " << elapsed << " s\n";
}

