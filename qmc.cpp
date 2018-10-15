#include <iostream>
#include "mcstep.hpp"

int main(){
    // Declare system
    qmc::MCStep engine;
    for (int i=0; i<19; i++){
#ifdef DEBUG_MODE
        std::cout << "Sweep " << i << std::endl;
#endif // DEBUG_MODE
        engine.DiagUpdate();
    }

    engine.LinkVertices();
#ifdef DEBUG_MODE
    std::cout << "Ok" << std::endl;
#endif // DEBUG_MODE
    return 0;
}
