#include "mcstep.hpp"
#include <iostream>

using namespace qmc;

//! MCStep constructor
/*!
 * This constructor currently uses a set of hard-coded values for simplicity.
 * \param lattice Pointer to the  lattice chosen by the user.
 */
MCStep::MCStep(Lattice* const lattice) : lattice(lattice),beta(20), expansion_factor(0.2),rdm_engine(std::mt19937(5)) {

    int op_str_length = 6; // Number of time slots
    measurement_frequency_ = 1000;
    LoopMethod = "wolff";
    sse_data = new SSEData(op_str_length,lattice,&rdm_engine);
    // Initialize number generators
    random_bond = std::uniform_int_distribution<unsigned int>(0,lattice->GetBondNumber()-1);
    random_prob = std::uniform_real_distribution<double>(0.0,1.0);
    random_leg = std::uniform_int_distribution<unsigned int>(0,4*op_str_length-1);

    observables_ = std::vector<Observable*>();
    output_manager_ = new OManager(1,this);
}

//! Class destructor
/*!
 * Empty for now. Will be used when dynamic allocation is employed in the future
 * of the project.
 */
MCStep::~MCStep(){
    // Delete random number generator
    delete sse_data;
    delete output_manager_;
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
 * \return Inverse temperature
 */
double MCStep::GetBeta() const { return beta; }

/*! \brief Adds an observable to the system
 *  Pushes back an observable. Nothing more.
 *  \param obsrv Pointer to the mentionned observable.
 */
void MCStep::AddObservable(Observable* obsrv) {
    observables_.push_back(obsrv);
}

const std::vector<Observable*>& MCStep::GetObservables() const {
    return observables_;
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
    output_manager_->PrintInitLine();
    for(unsigned int idx_sweep(0); idx_sweep<num_sweep; idx_sweep++) {
        LOG(2,"Sweep " << idx_sweep << std::endl);
        LOG(2,"=========================" << std::endl);
        DiagUpdate();
        LinkVertices();
        // OperatorLoop is implicitely called from LinkVertices.
        sse_data->Equilibrate(expansion_factor);
        random_leg = std::uniform_int_distribution<unsigned int>(0,4*sse_data->GetOperatorNumber(false)-1);
        UpdateObservables();
        if(idx_sweep%measurement_frequency_==0 && idx_sweep) {
            output_manager_->PrintObservables(idx_sweep);
            for(auto it=observables_.begin(); it!=observables_.end(); ++it)
                (*it)->ClearObservable();
        }
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
    for(size_t i=0; i<sse_data->GetOperatorNumber(false); i++){
        if(sse_data->GetOperator(i)==-1){ // No operator here
            unsigned int bond_num(random_bond(rdm_engine));
            std::vector<int> random_bond_spins(sse_data->GetBondSpins(bond_num));
            if((random_bond_spins[0] != random_bond_spins[1]) && AuthInsert()){
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
    bool all_links_empty(true);

    // Link operators from 0->L-1
    for(unsigned int i=0; i<sse_data->GetOperatorNumber(false); i++){
        if(sse_data->GetOperator(i) != -1){
            all_links_empty = false;
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

    this->OperatorLoop(pos_first,all_links_empty);
}

//! Selects the pre-configured operator loop to use.
/*!
 * Possible operator loops are, currently:
 * 1. Swendsen-Wang (qmc::MCStep::LoopMethod="swenwang")
 * 2. Wolff (qmc::MCStep::LoopMethod="wolff")
 * \param pos_first Vector containing the leg index where a spin is first used.
 * \param all_links_empty Indicates no operators are inserted.
 */
void MCStep::OperatorLoop(std::vector<int> const& pos_first, bool all_links_empty) {
    if(LoopMethod=="swenwang")
        OperatorLoopSwenWang(pos_first);
    else if (LoopMethod=="wolff")
        OperatorLoopWolff(pos_first,all_links_empty);
}

//! Browses operator loops and flips them (Swendsen-Wang Method)
/*!
 * This method follows the linked vertex list created through
 * qmc::MCStep::LinkVertices to map operator loops. Occasionally, the spins
 * along each loop will be flipped (which may in some case mean the associated
 * operator type changes). Single spins define a loop on their own and might
 * also be flipped with probability 1/2.
 * \param pos_first Vector containing the leg index where a spin is first used.
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
 * \param pos_first Vector containing the leg index where a spin is first used.
 * \param all_links_empty Indicates no operators are inserted.
 */
void MCStep::OperatorLoopWolff(std::vector<int> const& pos_first, bool all_links_empty){
    LOG(2,"Entering OperatorLoop" << std::endl << std::endl);

    unsigned int leg_counter(0); // Counts legs to exit following loops

    if(not all_links_empty){
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
    }

    // Second loop: Flip free spins
    for (unsigned int idx_spin(0); idx_spin < pos_first.size(); idx_spin++) {
        if (pos_first[idx_spin] == -1 && random_prob(rdm_engine) < 0.5)
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
 * \return Boolean allowing (or not) a insertion
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
 * \return Boolean allowing (or not) a removal
 */
bool MCStep::AuthRemove() {
    unsigned int num_bond(lattice->GetBondNumber());
    double threshold_remove(std::min(2.0 * (sse_data->GetOperatorNumber(false)-sse_data->GetOperatorNumber()+1.) / (beta*num_bond),1.0));

    if(random_prob(rdm_engine) < threshold_remove)
        return true;
    return false;
}

//! Gets the spin indexes of a leg and its neighbour
/*!
 * This function relies on a chain conversion from leg to operator, to bond and
 * finally to spin indexes. The first two ae implicitely defined in the
 * program's storage architecture, the last relies on Lattice::GetBondSpins.
 * \param idx_leg Index of the leg
 * \return The index of the spins attached to the operator using idx_leg
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
