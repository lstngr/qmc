#include "ssedata.hpp"
#include <iostream>

using namespace qmc;

//! SSEData constructor
/*!
 * This constructor will be called from qmc::MCStep::MCStep() when initializing
 * the simulation's storage.
 */
SSEData::SSEData(unsigned int time_slots, unsigned int spin_number, std::mt19937* rdm_engine)
    : num_time_slots(time_slots)
{

    // Prepare random spin generation
    std::uniform_int_distribution<unsigned int> random_spin(0,1);

    spins = std::vector<int>(spin_number,0);
    oprtr = std::vector<int>(num_time_slots,-1);
    bond_vertex = std::vector<std::vector<int>>(num_time_slots,std::vector<int>(4,-1));

    for(unsigned int idx_spin=0; idx_spin<spin_number; idx_spin++)
        spins[idx_spin] = 2*random_spin(*rdm_engine) - 1;

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
unsigned int SSEData::GetOperatorNumber(bool inserted){
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

//! Returns two spins belonging to a given bond
/*!
 * Returns a two-elements vector containing the spin values of the given bond
 * index.
 * \param num_bond Number of the bond. Should be < than the number of spins.
 */
std::vector<int> SSEData::GetBondSpins(unsigned int num_bond){
    std::vector<int> spins(2,0);
    spins[0] = this->spins[num_bond];
    spins[1] = this->spins[(num_bond+1)%this->spins.size()];

    return spins;
}

//! Flips the spins in a given bond
/*!
 * Looks up the spins belonging to a given bond and flips them. Useful when
 * encountering an off-diagonal operator during a MC sweep.
 * \param num_bond Index of the bond. Should be < than the number of spins.
 */
void SSEData::FlipSpins(unsigned int num_bond){
    spins[num_bond] *= -1;
    spins[(num_bond+1)%this->spins.size()] *= -1;
}

//! Prints operators for debbuging purposes
void SSEData::PrintOperators(){
    for(auto&& val_oprtr : oprtr)
        std::cout << val_oprtr << " ";
    std::cout << std::endl;
}

//! Prints spins for debbuging purposes
void SSEData::PrintSpins(){
    for(auto&& val_spin : spins)
        std::cout << val_spin << " ";
    std::cout << std::endl;
}


unsigned int SSEData::GetSpinNumber(){
    return spins.size();
}
