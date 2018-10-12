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
    //TODO Boundaries included?
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

//! Set Beta parameter
/*!
 * Sets the thermodynamic Beta parameter. Can be called at anytime, but
 * specifically by the class constructor.
 */
void MCStep::SetBeta(double beta){
    this->beta = beta;
}

//! Performs a diagonal update on the system
/*!
 * Inserts or removes diagonal operators from the system. Also updates spins if
 * an off-diagonal operator is encountered.
 */
void MCStep::DiagUpdate(){

    std::cout << "Before: " << std::endl;
    sse_data->PrintOperators();
    sse_data->PrintSpins();

    // Iterate over all possible operator positions
    for(unsigned int i=0; i<sse_data->GetOperatorNumber(false); i++){
        if(sse_data->GetOperator(i)==-1){ // No operator here
            unsigned int bond_num(random_bond(this->rdm_engine));
            std::vector<int> random_bond_spins(sse_data->GetBondSpins(bond_num));
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
        else{ // Off-Diagonal operator, update spins
            unsigned int bond = (sse_data->GetOperator(i)-1)/2;
            sse_data->FlipSpins(bond);
        }
    } // for over op_chain

    std::cout << "After: ";
    sse_data->PrintOperators();
    sse_data->PrintSpins();
    std::cout << std::endl;

} // MCStep::DiagUpdate


//! Diagonal insertion probability
/*!
 * Computes the probability of inserting a diagonal operator (from the current
 * state of the system).
 */
bool MCStep::AuthInsert() {
    unsigned int num_bond(sse_data->GetSpinNumber());
    double threshold_accept(std::min(beta*num_bond / (2.*(sse_data->GetOperatorNumber(false)-sse_data->GetOperatorNumber())),1.0));

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
    unsigned int num_bond(sse_data->GetSpinNumber());
    double threshold_remove(std::min(2.0 * (sse_data->GetOperatorNumber(false)-sse_data->GetOperatorNumber()+1) / (beta*num_bond),1.0));

    if(random_prob(rdm_engine) < threshold_remove)
        return true;
    return false;
}
