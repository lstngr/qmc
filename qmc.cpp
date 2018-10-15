#include <iostream>
#include "mcstep.hpp"

int main(){
    // Declare system
    qmc::MCStep engine;
    for (int i=0; i<19; i++){
        DEBUG_GEN("Sweep " << i << std::endl);
        engine.DiagUpdate();
    }

    engine.LinkVertices();
    DEBUG_GEN("Ok" << std::endl);
    return 0;
}
