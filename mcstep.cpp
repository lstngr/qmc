#include "mcstep.hpp"
#include <iostream>

using namespace qmc;

//! MCStep constructor
/*!
 * This constructor currently uses a set of hard-coded values for simplicity.
 */
MCStep::MCStep() : beta(10) {

    int spin_number = 6; // Numer of spin to simulate
    int op_str_length = 6; // Number of time slots
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

#ifdef DEBUG_GEN_BLOCK
    std::cerr << "Before: " << std::endl;
    sse_data->PrintOperators();
    sse_data->PrintSpins();
#endif //DEBUG_GEN_BLOCK

    // Iterate over all possible operator positions
    for(unsigned int i=0; i<sse_data->GetOperatorNumber(false); i++){
        if(sse_data->GetOperator(i)==-1){ // No operator here
            unsigned int bond_num(random_bond(this->rdm_engine));
            std::vector<int> random_bond_spins(sse_data->GetBondSpins(bond_num));
            if(random_bond_spins[0] != random_bond_spins[1] && AuthInsert()){
                // Try insertion
                DEBUG_GEN("Diagonal Insertion (bond " << bond_num << ", spins " << random_bond_spins[0] << " " << random_bond_spins[1] << ")" << std::endl);
                sse_data->SetOperator(i, 2 * bond_num);
            }
        }
        else if(sse_data->GetOperator(i)%2==0){ // There exists a diagonal operator
            // try remove
            if(AuthRemove()){
                DEBUG_GEN("Diagonal Removal (operator " << i << ")" << std::endl);
                sse_data->SetOperator(i,-1);
            }
        }
        else{ // Off-Diagonal operator, update spins
            unsigned int bond = sse_data->GetOperator(i)/2;
            sse_data->FlipSpins(bond);
        }
    } // for over op_chain

#ifdef DEBUG_GEN_BLOCK
    std::cerr << "After: ";
    sse_data->PrintOperators();
    sse_data->PrintSpins();
    std::cerr << std::endl;
#endif //DEBUG_GEN_BLOCK

} // MCStep::DiagUpdate


//! Updates the linked vertex list
/*!
 * Between operator interactions, each spin must be connected with it's
 * imaginary time partner. This function iterates on a "linked vertex" list in
 * order to create those references.
 */
void MCStep::LinkVertices(){

    std::vector<int> pos_last (sse_data->GetSpinNumber(),-1);
    std::vector<int> pos_first(sse_data->GetSpinNumber(),-1);

    // Link operators from 0->L-1
    for(unsigned int i=0; i<sse_data->GetOperatorNumber(false); i++){
        if(sse_data->GetOperator(i) != -1){
            int pos_oprtr(4*i);
            int idx_spin_1(sse_data->GetOperator(i)/2);
            int idx_spin_2(sse_data->GetOperator(i)/2+1);

            DEBUG_GEN("pos_oprtr=" << pos_oprtr << ", spins " << idx_spin_1 << " " << idx_spin_2 << std::endl);

            if(pos_last[idx_spin_1] != -1)
                sse_data->LinkVertex(pos_last[idx_spin_1],pos_oprtr);
            else
                pos_first[idx_spin_1] = pos_oprtr;

            if(pos_last[idx_spin_2] != -1)
                sse_data->LinkVertex(pos_last[idx_spin_2],pos_oprtr+1);
            else
                pos_first[idx_spin_2] = pos_oprtr + 1;

            pos_last[idx_spin_1] = pos_oprtr + 2;
            pos_last[idx_spin_2] = pos_oprtr + 3;
        }
    }

    // Link across last time boundary L-1->0
    for(unsigned int i=0; i<pos_first.size(); i++){
        if(pos_first[i] != -1)
            sse_data->LinkVertex(pos_first[i],pos_last[i]);
    }

#ifdef DEBUG_GEN_BLOCK
    sse_data->PrintLinks();
    std::cerr << std::endl;
#endif //DEBUG_GEN_BLOCK
}

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
