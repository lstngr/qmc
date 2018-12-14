#include "ssedata.hpp"
#include "lattice.hpp"
#include <iostream>
#include <iomanip>

using namespace qmc;

//! SSEData constructor
/*!
 * This constructor will be called from qmc::MCStep::MCStep() when initializing
 * the simulation's storage.
 * \param time_slots Initial number of time slots to create
 * \param lattice Pointer to the user's lattice, passed by qmc::MCStep::MCStep
 * \param rdm_engine Pointer to random number generator, shared with qmc::MCStep
 */
SSEData::SSEData(unsigned int time_slots, Lattice* const lattice, std::mt19937* const rdm_engine)
    : num_time_slots(time_slots), num_oprtr_type_(lattice->HowManyOperators()), pos_cfg(0), neg_cfg(0), lattice(lattice)
{
    // Prepare random spin generation
    spins = std::vector<Spin>(lattice->GetSpinNumber(),Spin(-1,0));
    oprtr = std::vector<int>(num_time_slots,-1);
    link_vertex = std::vector<int>(4*num_time_slots,-1);
    bond_vertex_= std::vector<std::vector<Spin>>(num_time_slots, std::vector<Spin>(4,Spin(-1,0)));

    for(unsigned int idx_spin=0; idx_spin<lattice->GetSpinNumber(); idx_spin++){
        // TODO - Add spins with multiple S?
        spins[idx_spin] = lattice->RandomSpin();
    }
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
 * \return The number of (inserted) operators in the system.
 */
unsigned int SSEData::GetOperatorNumber(bool inserted) const{
    if(not inserted)
        return oprtr.size();

    int empty_oprtr(0);
    for(unsigned int idx=0;idx<oprtr.size();idx++){
        if(oprtr[idx]==-1)
            empty_oprtr++;
    }

    return oprtr.size() - empty_oprtr;
}

//! Returns the specified spin value
//! \param idx_spin Index of the spin
//! \return Spin value at index idx_spin
Spin SSEData::GetSpin(unsigned int idx_spin) const{
    return spins[idx_spin];
}

//! Sets the specified spin to a given value
//! \param idx_spin Index of the spin
//! \param val_spin Value of the spin
void SSEData::SetSpin(unsigned int idx_spin, Spin spn){
    spins[idx_spin] = spn;
}

void SSEData::SetSpin(unsigned int idx_spin, int two_s, int two_sz){
    SetSpin(idx_spin, Spin(two_s,two_sz));
}

void SSEData::SetSpinVrtx(unsigned int idx_leg, Spin spn){
    bond_vertex_[idx_leg/4][idx_leg%4] = spn;
}

Spin SSEData::GetLegSpin(unsigned int leg_idx) const{
    return bond_vertex_[leg_idx/4][leg_idx%4];
}

std::vector<Spin> SSEData::GetBondVrtx(unsigned int oprtr_idx) const{
    return bond_vertex_[oprtr_idx];
}

void SSEData::RaiseSpin(unsigned int idx_leg, int amount){
    bond_vertex_[idx_leg/4][idx_leg%4].Raise(amount);
}

void SSEData::LowerSpin(unsigned int idx_leg, unsigned int amount){
    RaiseSpin(idx_leg,-amount);
}

//! Get the value of an operator (type and position).
/*!
 * Returns the s(p) = 2*b+a-1 value of the operator chain at a given
 * location (time slot).
 * \param idx_oprtr Time slot index where the operator lives.
 * \return Value of the operator at position idx_oprtr.
 */
int SSEData::GetOperator(unsigned int idx_oprtr) const{
    // Should implement index checks to avoid segmentation faults.
    return oprtr[idx_oprtr];
}

int SSEData::GetOperatorType(unsigned int idx_oprtr) const{
    return oprtr[idx_oprtr]%num_oprtr_type_;
}

int SSEData::GetOperatorBond(unsigned int idx_oprtr) const{
    return oprtr[idx_oprtr]/num_oprtr_type_;
}
//! Set the value of an operator at a given time slot
/*!
 * Sets the value of an operator at a given time slot index. If called without
 * setting value, the identity operator is inserted.
 * \param idx_oprtr The index (p) of the operator to set.
 * \param val_oprtr s(p)=2*b+a-1 value to set the given operator to.
 */
void SSEData::SetOperator(unsigned int idx_oprtr, int oprtr_type, int oprtr_bond, bool update_bond_vrtx){
    // Should implement index and value checks
    // Note: bool mainly for diagonal insertions/removals
    oprtr[idx_oprtr] = oprtr_type + num_oprtr_type_ * oprtr_bond;
    if(update_bond_vrtx){
        if(oprtr_type==-1)
            bond_vertex_[idx_oprtr] = std::vector<Spin>(4,Spin(-1,0));
        else if(lattice->IsDiagonal(oprtr_type)){
            unsigned int thebond(GetOperatorBond(idx_oprtr));
            // Okay while diag update since alpha being propagated
            std::vector<Spin> d_spins(this->GetBondSpins(thebond));
            d_spins.insert(std::end(d_spins), std::begin(d_spins), std::end(d_spins));
            bond_vertex_[idx_oprtr] = d_spins;
        }
    }
}

//! Switches the operator type at idx_oprtr
/*!
 * Carries out the conversion diagonal<->off-diagonal for the operator at the
 * given index.
 * \param idx_oprtr
 */
void SSEData::FlipOperator(unsigned int idx_oprtr) {
    SetOperator(idx_oprtr, GetOperator(idx_oprtr) ^ 1U << 0);
}

//! Returns a the link value at a given position
/*!
 * \param vertex Index of the link to fetch
 * \return Value of the link at vertex index.
 */
int SSEData::GetLink(unsigned int vertex) const {
    return link_vertex[vertex];
}

//! Set a value in the linked vertex storage
/*!
 * Modifies a link at a given position in the link_vertex storage. This function
 * is used when building operator loops and shouldn't be called otherwise.
 * \param vertex Index of the vertex to modify
 * \param val_vertex Value to assign
 */
void SSEData::SetLink(unsigned int vertex, int val_vertex) {
    link_vertex[vertex] = val_vertex;
}

//! Returns two spins belonging to a given bond
/*!
 * Returns a two-elements vector containing the spin values of the given bond
 * index.
 * \param num_bond Number of the bond. Should be < than the number of spins.
 * \return Vector containing the values of the spins linked by the bond with
 * index num_bond.
 */
std::vector<Spin> SSEData::GetBondSpins(unsigned int num_bond) const{
    std::vector<Spin> spins(2,Spin(-1,0));
    spins[0] = this->spins[num_bond];
    spins[1] = this->spins[(num_bond+1)%this->spins.size()];

    return spins;
}

void SSEData::UpdateSpins(unsigned int idx_oprtr){
    // TODO - Fix DIRTY stuff
    std::vector<Spin> curSpin(lattice->PropagateBond(*this,idx_oprtr));
    SetSpin(GetOperatorBond(idx_oprtr), curSpin[0]);
    SetSpin((GetOperatorBond(idx_oprtr)+1)%spins.size(), curSpin[1]);
}

//! Prints operators for debbuging purposes
void SSEData::PrintOperators() const{
    std::cout << "(" << GetOperatorNumber(false) << ") ";
    for(auto&& val_oprtr : oprtr)
        std::cout << val_oprtr << " ";
    std::cout << std::endl;
}

//! Prints spins for debbuging purposes
void SSEData::PrintSpins() const{
    for(auto&& val_spin : spins)
        std::cout << val_spin << " ";
    std::cout << std::endl;
}

//! Prints the vertex link table
void SSEData::PrintLinks() const{

    std::cout << "Linked Vertex Table" << std::endl;
    std::cout << std::left << std::setw(4) << "v " << std::left << std::setw(4) << "X(v)" << std::endl;

    for(unsigned int idx(0); idx<link_vertex.size(); idx++){
        std::cout << std::left << std::setw(4) << idx << std::left << std::setw(4) << link_vertex[idx] << std::endl;

    }

}

//! Creates a link between input vertices
/*!
 * Updates link_vertex so that the given indexes refer to each other.
 * \param vertex_0 One vertex index (4p+l)
 * \param vertex_1 Some other vertex index of the same form!
 */
void SSEData::LinkVertex(unsigned int vertex_0, unsigned int vertex_1){
    LOG(4,"Linking " << vertex_0 << " and " << vertex_1 << std::endl);
    link_vertex[vertex_0] = vertex_1;
    link_vertex[vertex_1] = vertex_0;
}

//! Resets the link_vertex to -1
void SSEData::ResetLinkVertex(){
    link_vertex = std::vector<int>(4*num_time_slots,-1);
    v_first_ = std::vector<int>(lattice->GetSpinNumber(),-1);
    v_last_  = std::vector<int>(lattice->GetSpinNumber(),-1);
}

int SSEData::GetVFirst(unsigned int idx_spin) const{
    return v_first_[idx_spin];
}

void SSEData::SetVFirst(unsigned int idx_spin, unsigned int idx_leg){
    v_first_[idx_spin] = idx_leg;
}

int SSEData::GetVLast(unsigned int idx_spin) const{
    return v_last_[idx_spin];
}

void SSEData::SetVLast(unsigned int idx_spin, unsigned int idx_leg){
    v_last_[idx_spin] = idx_leg;
}

//! Returns the next leg index for a SW&R procedure
/*!
 * The procedure of such an algorithm is to look up the neighbour spin (using a
 * bitflip to save CPU and memory), and then query its linked counterpart via
 * SSEData::GetLink().
 * \param idx_leg Index of the leg to start from.
 * \return Index of the leg after a SWR procedure.
 */
int SSEData::SwitchAndReverse(int idx_leg) const {
    // idx_leg XOR 1
    return idx_leg ^ 1 << 0;
}

int SSEData::ContinueStraight(int idx_leg) const {
    // idx_leg XOR 10
    return idx_leg ^ 1 << 1;
}

int SSEData::SwitchAndContinue(int idx_leg) const {
    return SwitchAndReverse(ContinueStraight(idx_leg));
}

void SSEData::UpdateSpinState(){
    for(unsigned int i=0;i<lattice->GetSpinNumber();i++){
        int first_leg(v_first_[i]);
        if(first_leg!=-1)
            spins[i] = bond_vertex_[first_leg/4][first_leg%4];
    }
}

//! Expands the operator time slots until the system equilibriates
/*!
 * Controls the expansion of the oprtr vector. The expansion factor is obtained
 * from MCStep.
 * \param expansion_factor Sets the minimum size of the lattice to L=(1+a)*n
 * where a is the discussed parameter.
 */
void SSEData::Equilibrate(double expansion_factor) {
    double thres_oprtr((1+expansion_factor) * GetOperatorNumber());
    if (num_time_slots < thres_oprtr) {
        // The system isn't large enough, expand it by adding unity operators.
        int delta_oprtr(std::ceil(thres_oprtr-num_time_slots));
        num_time_slots += delta_oprtr;
        std::vector<int> new_slots(delta_oprtr,-1);
        std::vector<std::vector<Spin>> new_vrtx(delta_oprtr,std::vector<Spin>(4,Spin(-1,0)));
        oprtr.insert(oprtr.end(),new_slots.begin(),new_slots.end());
        bond_vertex_.insert(bond_vertex_.end(),new_vrtx.begin(),new_vrtx.end());
    }
}

bool SSEData::CheckPositiveWeight(){
    int neg_count(0);
    for(unsigned int idx_oprtr(0); idx_oprtr<oprtr.size(); idx_oprtr++)
        neg_count += lattice->IsNegativeOperator(*this,4*idx_oprtr);
    if(neg_count%2){
        neg_cfg+=1;
        return false;
    }
    pos_cfg+=1;
    return true;
}

double SSEData::GetNegativeWeightRatio() {
    double ratio(1.*neg_cfg/(pos_cfg+neg_cfg));
    pos_cfg=0; neg_cfg=0;
    return ratio;
}
