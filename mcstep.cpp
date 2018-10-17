#include "mcstep.hpp"
#include <iostream>

using namespace qmc;

//! MCStep constructor
/*!
 * This constructor currently uses a set of hard-coded values for simplicity.
 */
MCStep::MCStep(Lattice const& lattice) : lattice(lattice),beta(10) {

    int spin_number = this->lattice.GetSpinNumber();
    int op_str_length = 100; // Number of time slots
    rdm_engine = std::mt19937(7); // Seed is 0
    sse_data = new SSEData(op_str_length,spin_number,&rdm_engine);
    // Initialize number generators
    random_bond = std::uniform_int_distribution<unsigned int>(0,lattice.GetBondNumber()-1);
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

    std::vector<int> pos_last (lattice.GetSpinNumber(),-1);
    std::vector<int> pos_first(lattice.GetSpinNumber(),-1);

    // Link operators from 0->L-1
    for(unsigned int i=0; i<sse_data->GetOperatorNumber(false); i++){
        if(sse_data->GetOperator(i) != -1){
            int pos_oprtr(4*i);
            std::vector<unsigned int> idx_spin(lattice.GetSpinsIndex(
                        sse_data->GetOperator(i)/2
                        ));

            DEBUG_GEN("pos_oprtr=" << pos_oprtr << ", spins " << idx_spin[0] << " " << idx_spin[1] << std::endl);

            if(pos_last[idx_spin[0]] != -1)
                sse_data->LinkVertex(pos_last[idx_spin[0]],pos_oprtr);
            else
                pos_first[idx_spin[0]] = pos_oprtr;

            if(pos_last[idx_spin[1]] != -1)
                sse_data->LinkVertex(pos_last[idx_spin[1]],pos_oprtr+1);
            else
                pos_first[idx_spin[1]] = pos_oprtr + 1;

            pos_last[idx_spin[0]] = pos_oprtr + 2;
            pos_last[idx_spin[1]] = pos_oprtr + 3;
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

    this->OperatorLoop(pos_first);
}

//! Browses operator loops and flips them
/*!
 * This method follows the linked vertex list created through
 * qmc::MCStep::LinkVertices to map operator loops. Occasionally, the spins
 * along each loop will be flipped (which may in some case mean the associated
 * operator type changes). Single spins define a loop on their own and might
 * also be flipped with probability 1/2.
 */
void MCStep::OperatorLoop(std::vector<int> const& pos_first){
    DEBUG_GEN("Entering OperatorLoop" << std::endl << std::endl);

    // First loop: Decides which oprtrs/spins to flip
    for (unsigned int idx_leg(0); idx_leg<4*sse_data->GetOperatorNumber(false); idx_leg+=2) {
        // Needs to be int for while exit criterion
        int idx_loop(idx_leg);
        DEBUG_GEN("Just set idx_loop=" << idx_loop << std::endl);
        if (sse_data->GetLink(idx_loop) >= 0) {
            // Determine if loop will be flipped:
            // -1 will not be, -2 will be,
            // Note these happen with equal probabilities (see random_prob)
            int val_loop(-1 - std::round(random_prob(rdm_engine)));
            DEBUG_GEN("Initial idx_loop=" << idx_loop << " with val_loop=" << val_loop << std::endl);
            do{
                DEBUG_GEN("Set val_loop in " << sse_data->GetLink(idx_loop) << "/" << idx_loop << std::endl);
                sse_data->SetLink(sse_data->GetLink(idx_loop),val_loop); // Sets 'previous' link to val_loop
                sse_data->SetLink(idx_loop,val_loop); // Sets own since link no longer needed
                sse_data->SetOperator(idx_loop/4, sse_data->GetOperator(idx_loop/4) ^ 1UL << 0); // Last bit flip to indicate diag<->off-diag
                idx_loop ^= 1UL << 0; // Flip last bit (ie looks for same side neighboor)
                DEBUG_GEN("Switch to neighbour idx_loop=" << idx_loop << std::endl);
                idx_loop = sse_data->GetLink(idx_loop);
                DEBUG_GEN("Goto next idx_loop=" << idx_loop << std::endl << std::endl);
            } while (idx_loop >= 0);
            // The loop exists when idx_loop points to a now modified location.
        }
    }

    // Second loop: Flip spins that are tagged either -1 or -2 using pos_first
    for (unsigned int idx_spin(0); idx_spin < pos_first.size(); idx_spin++) {
        int val_vertex(0);
        if (pos_first[idx_spin] != -1)
            val_vertex = sse_data->GetLink(pos_first[idx_spin]);
        else
            val_vertex = -1;
        if ((pos_first[idx_spin] == -1 && std::round(random_prob(rdm_engine))) || val_vertex == -2 ) {
            // If spin is marked as single and with a 50% chance,
            // OR if the attached vertex value mentions a flip,
            // Flip the spin!
            sse_data->SetSpin(idx_spin, -sse_data->GetSpin(idx_spin));
        }
    }

#ifdef DEBUG_GEN_BLOCK
    std::cerr << "After OperatorLoop: ";
    sse_data->PrintOperators();
    sse_data->PrintSpins();
    std::cerr << std::endl;
#endif //DEBUG_GEN_BLOCK

}

//! Diagonal insertion probability
/*!
 * Computes the probability of inserting a diagonal operator (from the current
 * state of the system).
 */
bool MCStep::AuthInsert() {
    unsigned int num_bond(lattice.GetBondNumber());
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
    unsigned int num_bond(lattice.GetSpinNumber());
    double threshold_remove(std::min(2.0 * (sse_data->GetOperatorNumber(false)-sse_data->GetOperatorNumber()+1) / (beta*num_bond),1.0));

    if(random_prob(rdm_engine) < threshold_remove)
        return true;
    return false;
}
