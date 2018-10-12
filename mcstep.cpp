#include "mcstep.hpp"
#include <iostream>

using namespace qmc;

//! MCStep constructor
/*!
 * This constructor currently uses a set of hard-coded values for simplicity.
 */
MCStep::MCStep() : beta(10) {

    int spin_number = 8; // Numer of spin to simulate
    int op_str_length = 16; // Number of time slots 
    rdm_engine = std::mt19937(1); // Seed is 0
    sse_data = new SSEData(op_str_length,spin_number,&rdm_engine);

    // Initialize number generators
    random_bond = std::uniform_int_distribution<unsigned int>(0,spin_number-1);
    random_prob = std::uniform_real_distribution<double>(0.0,1.0);
}

//! Class destructor
/*!
 * Empty for now. Will be used when dynamic allocation is employed in the future
 * of the project.
 */
MCStep::~MCStep(){
    // Close random number generator
    delete sse_data;
}

//! Performs a diagonal update on the system
/*!
 * Inserts or removes diagonal operators from the system. Also updates spins if
 * an off-diagonal operator is encountered.
 */
void MCStep::DiagUpdate(){

    std::cout << "Before: " << std::endl;
    this->PrintOperators();
    this->PrintSpins();

    // Iterate over all possible operator positions
    for(unsigned int i=0; i<sse_data->OperatorNumber(false); i++){
        if(sse_data->GetOperator(i)==-1){ // No operator here
            unsigned int bond_num(random_bond(this->rdm_engine));
            std::vector<int> random_bond_spins(this->BondToSpins(bond_num));
            if(random_bond_spins[0] != random_bond_spins[1] && AuthInsert()){
                // Try insertion
                std::cout << "Diagonal Insertion (bond " << bond_num << ", spins " << random_bond_spins[0] << " " << random_bond_spins[1] << ")" << std::endl;
                sse_data->SetOperator(i, 2 * bond_num);
            }
        }
        else if(sse_data->GetOperator(i)%2==0){ // There exists a diagonal operator
            // try remove
            if(AuthRemove()){
                std::cout << "Diagonal Removal (operator " << i << ")" << std::endl;
                sse_data->SetOperator(i,-1);
            }
        }
        else if(sse_data->GetOperator(i)!=-1){ // Off-Diagonal operator, update spins
            unsigned int bond = (sse_data->GetOperator(i)-1)/2;
            this->FlipSpins(bond);
        }
    } // for over op_chain

    std::cout << "After: ";
    this->PrintOperators();
    this->PrintSpins();
    std::cout << std::endl;

} // MCStep::DiagUpdate


//! Takes a bond index, returns the two corresponding spin indexes
std::vector<int> MCStep::BondToSpins(unsigned int num_bond){
    std::vector<int> spins(2,0);
    spins[0] = this->spins[num_bond];
    spins[1] = this->spins[(num_bond+1)%this->spins.size()]; // Works for a 1D circular chain
    return spins;
} // MCStep::BondToSpins


//! Flips the spins belonging to a given bond
void MCStep::FlipSpins(unsigned int num_bond){
    this->spins[num_bond] = - this->spins[num_bond];
    this->spins[(num_bond+1)%this->spins.size()] = - this->spins[(num_bond+1)%this->spins.size()];
} // MCStep::FlipSpins


//! Displays the operator chain for debugging purposes
void MCStep::PrintOperators(){
    for (auto op : this->op_chain)
        std::cout << op << " ";
    std::cout << std::endl;
}

//! Displays the operator chain for debugging purposes
void MCStep::PrintSpins(){
    for (auto sp : this->spins)
        std::cout << sp << " ";
    std::cout << std::endl;
}

//! Diagonal insertion probability
/*!
 * Computes the probability of inserting a diagonal operator (from the current
 * state of the system).
 */
bool MCStep::AuthInsert() {
    unsigned int num_bond(spins.size());
    double threshold_accept(std::min(beta*num_bond / (2.*(op_str_length-num_cur_op)),1.0));

    if(random_prob(rdm_engine) < threshold_accept)
        return true;
    return false;
}

//! Diagonal removal probability
/*!
 * Computes the probability of removing a diagonal operator (from the current
 * state of the system).
 */
bool MCStep::AuthRemove() {
    unsigned int num_bond(spins.size());
    double threshold_remove(std::min(2.0 * (op_str_length-num_cur_op+1) / (beta*num_bond),1.0));

    if(random_prob(rdm_engine) < threshold_remove)
        return true;
    return false;
}
