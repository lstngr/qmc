#include "mcstep.hpp"
#include <iostream>

using namespace qmc;

//! MCStep constructor
/*!
 * This constructor currently uses a set of hard-coded values for simplicity.
 */
MCStep::MCStep(Lattice* const lattice) : lattice(lattice),beta(100), expansion_factor(0.4),rdm_engine(std::mt19937(7)) {

    int op_str_length = 6; // Number of time slots
    LoopMethod = "wolff";
    sse_data = new SSEData(op_str_length,lattice,&rdm_engine);
    // Initialize number generators
    random_bond = std::uniform_int_distribution<unsigned int>(0,lattice->GetBondNumber()-1);
    random_prob = std::uniform_real_distribution<double>(0.0,1.0);
    random_leg = std::uniform_int_distribution<unsigned int>(0,4*op_str_length-1);

    observables_ = std::vector<Observable*>();
}

//! Class destructor
/*!
 * Empty for now. Will be used when dynamic allocation is employed in the future
 * of the project.
 */
MCStep::~MCStep(){
    // Delete random number generator
    delete sse_data;
}

//! Set Beta parameter
/*!
 * Sets the thermodynamic Beta parameter. Can be called at anytime, but
 * specifically by the class constructor.
 * \param beta The value of the thermodynamic parameter.
 */
void MCStep::SetBeta(double beta){
    this->beta = beta;
}

//! Get Beta parameter
/*!
 * Returns the set thermodynamic beta.
 */
double MCStep::GetBeta() const { return beta; }

/*! \brief Adds an observable to the system
 *  Pushes back an observable. Nothing more.
 *  \param obsrv Pointer to the mentionned observable.
 */
void MCStep::AddObservable(Observable* obsrv) {
    observables_.push_back(obsrv);
}

//! Run a full Monte-Carlo Sweep
/*!
 * This function runs a full quantum Monte-Carlo sweep on the system, that is
 * - A diagonal update (qmc::MCStep::DiagUpdate)
 * - Vertex legs linkage (qmc::MCStep::LinkVertices)
 * - Building operators loops and flipping them (qmc::MCStep::OperatorLoop)
 *
 * Due to the interdependence of these methods on each other, they are flagged
 * as private and called sequentially by the present function.
 *
 * \param num_sweep How many sweep should be carried out
 */
void MCStep::Sweep(unsigned int num_sweep){
    for(unsigned int idx_sweep(0); idx_sweep<num_sweep; idx_sweep++) {
        LOG(2,"Sweep " << idx_sweep << std::endl);
        LOG(2,"=========================" << std::endl);
        DiagUpdate();
        LinkVertices();
        // OperatorLoop is implicitely called from LinkVertices.
        sse_data->Equilibrate(expansion_factor);
        random_leg = std::uniform_int_distribution<unsigned int>(0,4*sse_data->GetOperatorNumber(false)-1);
        UpdateObservables();
    }
}

//! Performs a diagonal update on the system
/*!
 * Inserts or removes diagonal operators from the system. Also updates spins if
 * an off-diagonal operator is encountered.
 */
void MCStep::DiagUpdate(){

#if LOG_LEVEL > 1
    std::cerr << "Before: " << std::endl;
    sse_data->PrintOperators();
    sse_data->PrintSpins();
#endif // LOG_LEVEL > 1

    // Iterate over all possible operator positions
    for(unsigned int i=0; i<sse_data->GetOperatorNumber(false); i++){
        if(sse_data->GetOperator(i)==-1){ // No operator here
            unsigned int bond_num(random_bond(rdm_engine));
            std::vector<int> random_bond_spins(sse_data->GetBondSpins(bond_num));
            if(random_bond_spins[0] != random_bond_spins[1] && AuthInsert()){
                // Try insertion
                LOG(3,"Diagonal Insertion (bond " << bond_num << ", spins " << random_bond_spins[0] << " " << random_bond_spins[1] << ")" << std::endl);
                sse_data->SetOperator(i, 2 * bond_num);
            }
        }
        else if(sse_data->GetOperator(i)%2==0){ // There exists a diagonal operator
            // try remove
            if(AuthRemove()){
                LOG(3,"Diagonal Removal (operator " << i << ")" << std::endl);
                sse_data->SetOperator(i,-1);
            }
        }
        else{ // Off-Diagonal operator, update spins
            unsigned int bond = sse_data->GetOperator(i)/2;
            sse_data->FlipSpins(bond);
        }
    } // for over op_chain

#if LOG_LEVEL > 1
    std::cerr << "After: " << std::endl;
    sse_data->PrintOperators();
    sse_data->PrintSpins();
    std::cerr << std::endl;
#endif // LOG_LEVEL > 1

} // MCStep::DiagUpdate


//! Updates the linked vertex list
/*!
 * Between operator interactions, each spin must be connected with it's
 * imaginary time partner. This function iterates on a "linked vertex" list in
 * order to create those references.
 */
void MCStep::LinkVertices(){

    sse_data->ResetLinkVertex();
    std::vector<int> pos_last (lattice->GetSpinNumber(),-1);
    std::vector<int> pos_first(lattice->GetSpinNumber(),-1);

    // Link operators from 0->L-1
    for(unsigned int i=0; i<sse_data->GetOperatorNumber(false); i++){
        if(sse_data->GetOperator(i) != -1){
            int pos_oprtr(4*i);
            std::vector<unsigned int> idx_spin(lattice->GetSpinsIndex(
                        sse_data->GetOperator(i)/2
                        ));

            LOG(4,"pos_oprtr=" << pos_oprtr << ", spins " << idx_spin[0] << " " << idx_spin[1] << std::endl);

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

#if LOG_LEVEL > 2
    sse_data->PrintLinks();
    std::cerr << std::endl;
#endif // LOG_LEVEL > 1

    this->OperatorLoop(pos_first);
}

//! Selects the pre-configured operator loop to use.
/*!
 * Possible operator loops are, currently:
 * 1) Swendsen-Wang (qmc::MCStep::LoopMethod="swenwang")
 * 2) Wolff (qmc::MCStep::LoopMethod="wolff")
 */
void MCStep::OperatorLoop(std::vector<int> const& pos_first) {
    if(LoopMethod=="swenwang")
        OperatorLoopSwenWang(pos_first);
    else if (LoopMethod=="wolff")
        OperatorLoopWolff(pos_first);
}

//! Browses operator loops and flips them (Swendsen-Wang Method)
/*!
 * This method follows the linked vertex list created through
 * qmc::MCStep::LinkVertices to map operator loops. Occasionally, the spins
 * along each loop will be flipped (which may in some case mean the associated
 * operator type changes). Single spins define a loop on their own and might
 * also be flipped with probability 1/2.
 */
void MCStep::OperatorLoopSwenWang(std::vector<int> const& pos_first){
    LOG(2,"Entering OperatorLoop" << std::endl << std::endl);

    // First loop: Decides which oprtrs/spins to flip
    for (unsigned int idx_leg(0); idx_leg<4*sse_data->GetOperatorNumber(false); idx_leg+=2) {
        // Needs to be int for while exit criterion
        int idx_loop(idx_leg);
        LOG(4,"Just set idx_loop=" << idx_loop << std::endl);
        if (sse_data->GetLink(idx_loop) >= 0) {
            // Determine if loop will be flipped:
            // -1 will not be, -2 will be,
            // Note these happen with equal probabilities (see random_prob)
            int val_loop(-1 - std::round(random_prob(rdm_engine)));
            LOG(3,"Initial idx_loop=" << idx_loop << " with val_loop=" << val_loop << std::endl);
            do{
                LOG(3,"Set val_loop in " << sse_data->GetLink(idx_loop) << "/" << idx_loop << std::endl);
                sse_data->SetLink(sse_data->GetLink(idx_loop),val_loop); // Sets 'previous' link to val_loop
                sse_data->SetLink(idx_loop,val_loop); // Sets own since link no longer needed
                if (val_loop==-2) {
                    LOG(3,"Will set op:" << sse_data->GetOperator(idx_loop/4) << "->" << (sse_data->GetOperator(idx_loop/4)^1 << 0) << std::endl);
                    sse_data->FlipOperator(idx_loop/4); // Last bit flip to indicate diag<->off-diag
                }
                LOG(3,"Switch to neighbour idx_loop=" << (idx_loop ^ 1 << 0) << std::endl);
                idx_loop = sse_data->SwitchAndReverse(idx_loop);
                LOG(3,"Goto next idx_loop=" << idx_loop << std::endl << std::endl);
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
            sse_data->FlipSpin(idx_spin);
        }
    }

#if LOG_LEVEL > 1
    std::cerr << "After OperatorLoop: ";
    sse_data->PrintOperators();
    sse_data->PrintSpins();
    std::cerr << std::endl;
#endif // LOG_LEVEL > 1

}
//! Browses operator loops and flips them (Wolff Method)
/*!
 * This method follows the linked vertex list created through
 * qmc::MCStep::LinkVertices to map operator loops. Occasionally, the spins
 * along each loop will be flipped (which may in some case mean the associated
 * operator type changes). Single spins define a loop on their own and might
 * also be flipped with probability 1/2.
 */
void MCStep::OperatorLoopWolff(std::vector<int> const& pos_first){
    LOG(2,"Entering OperatorLoop" << std::endl << std::endl);

    unsigned int leg_counter(0); // Counts legs to exit following loops

    // First loop: Flip loops
    do {
        // Needs to be int for while exit criterion
        int idx_loop(random_leg(rdm_engine));
        int idx_loop_start(idx_loop);
        LOG(3,"Start set idx_loop=" << idx_loop << std::endl);
        if (sse_data->GetLink(idx_loop) >= 0) {
            do{
                sse_data->FlipOperator(idx_loop/4);
                LOG(3,"Set oprtr " << (idx_loop/4) << " from " << (sse_data->GetOperator(idx_loop/4) ^ 1 << 0) << " to " << sse_data->GetOperator(idx_loop/4) << std::endl);
                LOG(3,"Switch to neighbour " << (idx_loop ^ 1 << 0) << " and then link to " << sse_data->GetLink(idx_loop ^ 1 << 0) << std::endl);
                idx_loop = sse_data->SwitchAndReverse(idx_loop);

                // Check if flipped spins belong to pos_first
                std::vector<unsigned int> idx_spins(LegToSpins(idx_loop));
                LOG(4,"Got bond spins " << idx_spins[0] << " " << idx_spins[1] << std::endl);
                if (pos_first[idx_spins[0]] == idx_loop || pos_first[idx_spins[0]] == (idx_loop ^ 1 << 0)) {
                    sse_data->FlipSpin(idx_spins[0]);
                    LOG(4,"Flip spin (idx " << idx_spins[0] << ")." << std::endl);
                }
                if (pos_first[idx_spins[1]] == idx_loop || pos_first[idx_spins[1]] == (idx_loop ^ 1 << 0)) {
                    LOG(4,"Flip spin (idx " << idx_spins[1] << ")." << std::endl);
                    sse_data->FlipSpin(idx_spins[1]);
                }

                // Update leg counter
                leg_counter++;
            } while (idx_loop != idx_loop_start);
            // The loop exists when idx_loop points to the start location.
        }
    } while(leg_counter < 4*sse_data->GetOperatorNumber(false));
    // Stop when enough loops have been performed

    // Second loop: Flip free spins
    for (unsigned int idx_spin(0); idx_spin < pos_first.size(); idx_spin++) {
        if (pos_first[idx_spin] == -1)
            sse_data->FlipSpin(idx_spin);
    }

#if LOG_LEVEL > 1
    std::cerr << "After OperatorLoop: ";
    sse_data->PrintOperators();
    sse_data->PrintSpins();
    std::cerr << std::endl;
#endif // LOG_LEVEL > 1

}

//! Diagonal insertion probability
/*!
 * Computes the probability of inserting a diagonal operator (from the current
 * state of the system).
 */
bool MCStep::AuthInsert() {
    unsigned int num_bond(lattice->GetBondNumber());
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
    unsigned int num_bond(lattice->GetSpinNumber());
    double threshold_remove(std::min(2.0 * (sse_data->GetOperatorNumber(false)-sse_data->GetOperatorNumber()+1) / (beta*num_bond),1.0));

    if(random_prob(rdm_engine) < threshold_remove)
        return true;
    return false;
}

//! Gets the spin indexes of a leg and its neighbour
/*!
 * This function relies on a chain conversion from leg to operator, to bond and
 * finally to spin indexes. The first two ae implicitely defined in the
 * program's storage architecture, the last relies on Lattice::GetBondSpins.
 */
std::vector<unsigned int> MCStep::LegToSpins(unsigned int idx_leg) const {
    int idx_oprtr(idx_leg/4);
    int idx_bond(sse_data->GetOperator(idx_oprtr)/2);
    std::vector<unsigned int> idx_spins(lattice->GetSpinsIndex(idx_bond));
    return idx_spins;
}

/*! \brief Updates observable states
 * Calls every observable's update function, so they store relevant informations
 * for their final value.
 */
void MCStep::UpdateObservables() const {
    for(auto obsrv : observables_) {
        obsrv->UpdateObservable(sse_data);
    }
}
