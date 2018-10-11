#include <iostream>
#include "mcstep.hpp"

int main(){
    // Declare system
    qmc::MCStep engine;
    for (int i=0; i<19; i++){
        std::cout << "Sweep " << i << std::endl;
        engine.DiagUpdate();
    }
    std::cout << "Ok" << std::endl;
    return 0;
}
