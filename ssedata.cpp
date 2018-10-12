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


//! Returns the number of operators in system.
/*!
 * Returns the number of operators in the system. Can return either the number
 * of currently inserted operators, or the maximum number of operators allowed.
 * \param inserted If true, the maximum number of operators is returned.
 */
int SSEData::OperatorNumber(bool inserted){
    if(not inserted)
        return oprtr.size();

    int empty_oprtr(0);
    for(unsigned int idx=0;idx<oprtr.size();idx++){
        if(oprtr[idx]==-1)
            empty_oprtr++;
    }

    return oprtr.size() - empty_oprtr;
}

//! Get the value of an operator (type and position).
/*!
 * Returns the s(p) = 2*b+a-1 value of the operator chain at a given
 * location (time slot).
 * \param idx_oprtr Time slot index where the operator lives.
 */
int SSEData::GetOperator(unsigned int idx_oprtr){
    // Should implement index checks to avoid segmentation faults.
    return oprtr[idx_oprtr];
}

//! Set the value of an operator at a given time slot
/*!
 * Sets the value of an operator at a given time slot index. If called without
 * setting value, the identity operator is inserted.
 * \param val_oprtr s(p)=2*b+a-1 value to set the given operator to.
 */
void SSEData::SetOperator(unsigned int idx_oprtr, int val_oprtr){
    // Should implement index and value checks
    oprtr[idx_oprtr] = val_oprtr;
}
