#include "ssedata.hpp"

using namespace qmc;

//! SSEData constructor
/*! 
 * This constructor will be called from qmc::MCStep::MCStep() when initializing
 * the simulation's storage.
 */
SSEData::SSEData(): num_time_slots(10){
    
    spins = std::vector<int>(50,0);
    oprtr = std::vector<int>(num_time_slots,0);
    bond_vertex = std::vector<std::vector<int>>(num_time_slots,std::vector<int>(4,0));

}

//! SSEData destructor
/*!
 * Might be useful later. No purpose at the moment of writing.
 */
SSEData::~SSEData(){}
