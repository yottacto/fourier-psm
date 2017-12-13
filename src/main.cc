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

    auto toml = cpptoml::parse_file("config.toml");
    icesp::configuration::config config{toml};
    // base number of cores per direction
    int bncd{config.cases.begin()->second.bncd};
    // base number of points per core per direction
    int bnpcd{config.cases.begin()->second.bnpcd};
    int iteration{config.cases.begin()->second.iteration};

    fpsm::solver s(bncd, bnpcd);

    s.init(fpsm::default_f);
    s.iterate(iteration);
    s.print();

    MPI::Finalize();
}

