#include <iostream>
#include "mcstep.hpp"

int main(){
    // Declare system
    qmc::Lattice1D lattice(6);
    qmc::MCStep engine(&lattice);
    engine.Sweep(3);

    DEBUG_GEN("Ok" << std::endl);
    return 0;
}
