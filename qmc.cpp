#include <iostream>
#include "mcstep.hpp"

int main(){
    // Declare system
    qmc::Lattice1D lattice(8);
    qmc::MCStep engine(&lattice);
    engine.Sweep(60);

    LOG(3,"Ok" << std::endl);
    return 0;
}
