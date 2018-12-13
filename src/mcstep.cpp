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
    sse_data = new SSEData(op_str_length,lattice,&rdm_engine);
    // Initialize number generators
    random_bond = std::uniform_int_distribution<unsigned int>(0,lattice->GetBondNumber()-1);
    random_prob = std::uniform_real_distribution<double>(0.0,1.0);
    random_leg = std::uniform_int_distribution<unsigned int>(0,4*op_str_length-1);

    observables_ = std::vector<Observable*>();
    output_manager_ = new OManager(1,this);
}

/*! Builds MCStep using a configuration file.
 * It still relies on the user for Lattice and Observable initialization!
 * \param lattice Pointer to a user-initialized lattice
 * \param sim_conf_ Configuration parameters
 */
MCStep::MCStep(Lattice* const lattice, const Configuration& sim_conf_)
    :lattice(lattice), rdm_engine(std::mt19937(5))
{
    int op_str_length(6);
    beta = sim_conf_.beta;
    expansion_factor = sim_conf_.expansion_factor;
    measurement_frequency_ = sim_conf_.measurement_frequency;
    sse_data = new SSEData(op_str_length,lattice,&rdm_engine);
    // Initialize number generators
    random_bond = std::uniform_int_distribution<unsigned int>(0,lattice->GetBondNumber()-1);
    random_prob = std::uniform_real_distribution<double>(0.0,1.0);
    random_leg = std::uniform_int_distribution<unsigned int>(0,4*op_str_length-1);
    observables_ = std::vector<Observable*>();
    output_manager_ = new OManager(1,this,sim_conf_.output_dest);
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

//! Provides (constant) access to observables_'s data.
//! \return Constant reference to system observables
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
        sse_data->UpdateSpinState();
        UpdateObservables();
        if(idx_sweep%measurement_frequency_==1 && idx_sweep) {
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
        unsigned int bond_num(random_bond(rdm_engine));
        if(sse_data->GetOperator(i)==-1){ // No operator here
            std::vector<Spin> random_bond_spins(sse_data->GetBondSpins(bond_num));
            int diag_oprtr_type(lattice->FindOperatorType(random_bond_spins));
            if(diag_oprtr_type!=-1 && AuthInsert(lattice->FindOperatorWeight(diag_oprtr_type))){
                // Try insertion
                sse_data->SetOperator(i, diag_oprtr_type, bond_num,true);
            }
        }
        else if(lattice->IsDiagonal(sse_data->GetOperatorType(i))){
            if(AuthRemove(lattice->FindOperatorWeight(sse_data->GetOperatorType(i)))){
                sse_data->SetOperator(i,-1,0,true);
            }
        }
        else{ // Off-Diagonal operator, update spins
            // Could call every iteration, but so prevents
            // unnecessary computations.
            sse_data->UpdateSpins(i);
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
    bool all_links_empty(true);

    // Link operators from 0->L-1
    for(unsigned int i=0; i<sse_data->GetOperatorNumber(false); i++){
        if(sse_data->GetOperator(i) != -1){
            all_links_empty = false;
            int pos_oprtr(4*i);
            std::vector<unsigned int> idx_spin(lattice->GetSpinsIndex(
                        sse_data->GetOperatorBond(i)
                        ));

            LOG(4,"pos_oprtr=" << pos_oprtr << ", spins " << idx_spin[0] << " " << idx_spin[1] << std::endl);

            if(sse_data->GetVLast(idx_spin[0]) != -1)
                sse_data->LinkVertex(sse_data->GetVLast(idx_spin[0]),pos_oprtr);
            else
                sse_data->SetVFirst(idx_spin[0], pos_oprtr);

            if(sse_data->GetVLast(idx_spin[1]) != -1)
                sse_data->LinkVertex(sse_data->GetVLast(idx_spin[1]),pos_oprtr+1);
            else
                sse_data->SetVFirst(idx_spin[1], pos_oprtr + 1);

            sse_data->SetVLast(idx_spin[0], pos_oprtr + 2);
            sse_data->SetVLast(idx_spin[1], pos_oprtr + 3);
        }
    }

    // Link across last time boundary L-1->0
    for(unsigned int i=0; i<lattice->GetSpinNumber(); i++){
        if(sse_data->GetVFirst(i) != -1)
            sse_data->LinkVertex(sse_data->GetVFirst(i),sse_data->GetVLast(i));
    }

#if LOG_LEVEL > 2
    sse_data->PrintLinks();
    std::cerr << std::endl;
#endif // LOG_LEVEL > 1

    this->OperatorLoop(all_links_empty);
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
void MCStep::OperatorLoop(bool all_links_empty){
    LOG(2,"Entering OperatorLoop" << std::endl << std::endl);

    unsigned int leg_counter(0); // Counts legs to exit following loops

    if(not all_links_empty){
        // First loop: Flip loops
        do {
            // Needs to be int for while exit criterion
            int idx_loop(random_leg(rdm_engine));
            int idx_loop_start(idx_loop);

            if (lattice->IsLoopStart(*sse_data,idx_loop)) {
                Spin loop_type(lattice->GetLoopType(*sse_data,idx_loop));
                do{
                    // Update bond vertex
                    idx_loop = sse_data->GetLink(idx_loop);
                    // Get Moves
                    idx_loop = lattice->ChooseMove(*sse_data, idx_loop, loop_type);
                    unsigned int thebond(sse_data->GetOperatorBond(idx_loop/4));
                    sse_data->SetOperator(idx_loop/4,
                            lattice->FindOperatorType(*sse_data, idx_loop),
                            thebond);
                    // TODO -- Fix bond determination

                    // Update leg counter
                    leg_counter++;
                } while (idx_loop != idx_loop_start);
                // The loop exists when idx_loop points to the start location.
            }
        } while(leg_counter < 4*sse_data->GetOperatorNumber(false));
        // Stop when enough loops have been performed
    }

    // Second loop: Flip free spins
    for (unsigned int idx_spin(0); idx_spin < lattice->GetSpinNumber(); idx_spin++) {
        if (sse_data->GetVFirst(idx_spin) == -1 && random_prob(rdm_engine) < 0.5){}// Choose random spin value
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
bool MCStep::AuthInsert(double weight) {
    unsigned int num_bond(lattice->GetBondNumber());
    double threshold_accept(std::min(beta*num_bond*weight / (sse_data->GetOperatorNumber(false)-sse_data->GetOperatorNumber()),1.0));

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
bool MCStep::AuthRemove(double weight) {
    unsigned int num_bond(lattice->GetBondNumber());
    double threshold_remove(std::min((sse_data->GetOperatorNumber(false)-sse_data->GetOperatorNumber()+1.) / (beta*num_bond*weight),1.0));

    if(random_prob(rdm_engine) < threshold_remove)
        return true;
    return false;
}

/*! \brief Updates observable states
 * Calls every observable's update function, so they store relevant informations
 * for their final value.
 */
void MCStep::UpdateObservables() const {
    for(auto obsrv : observables_)
        obsrv->UpdateObservable(sse_data);
}
