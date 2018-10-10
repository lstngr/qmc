#include <iostream>
#include "qmcstep.hpp"

int main(){
    // Declare system
    qmc::QMCStep engine;
    for (int i=0; i<19; i++){
        std::cout << "Sweep " << i << std::endl;
        engine.DiagUpdate();
    }
    std::cout << "Ok" << std::endl;
    return 0;
}
