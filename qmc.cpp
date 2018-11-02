#include <iostream>
#include "mcstep.hpp"

int main(){
    // Declare system
    qmc::Lattice1D lattice(16);
    qmc::IManager read_input;
    qmc::MCStep engine(&lattice, read_input.GetConfig());
    qmc::ObservableHamilton energy(engine.GetBeta());
    engine.AddObservable(&energy);
    engine.Sweep(read_input.GetConfig().num_sweep);

    std::cout << "Energy is: " << energy.Measure() << std::endl;
    LOG(3,"Ok" << std::endl);
    return 0;
}
