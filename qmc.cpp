#include <iostream>
#include "mcstep.hpp"

int main(int argc, char** argv){
    std::string config_file;
    if(argc==1)
        config_file=argv[0];
    // Declare system
    qmc::Lattice1D lattice(16);
    qmc::IManager read_input(config_file);
    qmc::MCStep engine(&lattice, read_input.GetConfig());
    qmc::ObservableHamilton energy(engine.GetBeta());
    engine.AddObservable(&energy);
    engine.Sweep(read_input.GetConfig().num_sweep);

    std::cout << "Energy is: " << energy.Measure() << std::endl;
    LOG(3,"Ok" << std::endl);
    return 0;
}
