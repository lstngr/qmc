#include <iostream>
#include "mcstep.hpp"

int main(){
    // Declare system
    qmc::Lattice1D lattice(8);
    qmc::MCStep engine(&lattice);
    qmc::ObservableHamilton energy(engine.GetBeta());
    engine.AddObservable(&energy);
    engine.Sweep(7500);

    std::cout << "Energy is: " << energy.Measure() << std::endl;
    LOG(3,"Ok" << std::endl);
    return 0;
}
