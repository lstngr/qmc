#include "ssedata.hpp"
#include <iostream>
#include <iomanip>

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
    link_vertex = std::vector<int>(4*num_time_slots,-1);

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
int SSEData::GetSpin(unsigned int idx_spin) const{
    return spins[idx_spin];
}

//! Sets the specified spin to a given value
//! \param idx_spin Index of the spin
//! \param val_spin Value of the spin
void SSEData::SetSpin(unsigned int idx_spin, int val_spin){
    spins[idx_spin] = val_spin;
}

//! Get the value of an operator (type and position).
/*!
 * Returns the s(p) = 2*b+a-1 value of the operator chain at a given
 * location (time slot).
 * \param idx_oprtr Time slot index where the operator lives.
 */
int SSEData::GetOperator(unsigned int idx_oprtr) const{
    // Should implement index checks to avoid segmentation faults.
    return oprtr[idx_oprtr];
}

//! Set the value of an operator at a given time slot
/*!
 * Sets the value of an operator at a given time slot index. If called without
 * setting value, the identity operator is inserted.
 * \param idx_oprtr The index (p) of the operator to set.
 * \param val_oprtr s(p)=2*b+a-1 value to set the given operator to.
 */
void SSEData::SetOperator(unsigned int idx_oprtr, int val_oprtr){
    // Should implement index and value checks
    oprtr[idx_oprtr] = val_oprtr;
}

//! Returns a the link value at a given position
/*!
 * \param vertex Index of the link to fetch
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

//! Returns the next element in an operator loop
/*!
 * Returns the index of the next element in a operator loop. Behavior unknown
 * !!!
 * \param vertex Index of the current vertex.
 */
int SSEData::GetLoopNeighbour(unsigned int vertex){
    return link_vertex[vertex];
}

//! Returns two spins belonging to a given bond
/*!
 * Returns a two-elements vector containing the spin values of the given bond
 * index.
 * \param num_bond Number of the bond. Should be < than the number of spins.
 */
std::vector<int> SSEData::GetBondSpins(unsigned int num_bond) const{
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
void SSEData::PrintOperators() const{
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

//! Returns current number of spins
unsigned int SSEData::GetSpinNumber() const{
    return spins.size();
}

//! Creates a link between input vertices
/*!
 * Updates link_vertex so that the given indexes refer to each other.
 * \param vertex_0 One vertex index (4p+l)
 * \param vertex_1 Some other vertex index of the same form!
 */
void SSEData::LinkVertex(unsigned int vertex_0, unsigned int vertex_1){
    DEBUG_GEN("Linking " << vertex_0 << " and " << vertex_1 << std::endl);
    link_vertex[vertex_0] = vertex_1;
    link_vertex[vertex_1] = vertex_0;
}

//! Resets the link_vertex to -1
void SSEData::ResetLinkVertex(){
    link_vertex = std::vector<int>(4*num_time_slots,-1);
}
