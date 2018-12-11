#include "lattice.hpp"
#include "ssedata.hpp"

using namespace qmc;

////////////////////////////
// LATTICE /////////////////
////////////////////////////

//! Lattice class constructor
Lattice::Lattice(){

}


////////////////////////////
// LATTICE1D ///////////////
////////////////////////////

//! Lattice1D class constructor
/*!
 * Instantiates a linear structure containing spins. This structure is cyclic
 * (hence, the number of spins is also the number of bonds).
 * \param length Length of the lattice, ie. the site number.
 */
Lattice1D::Lattice1D(unsigned int length): length(length){}

//! Returns spin indexes of a specific bond
/*!
 * Returns the two indexes of each spin belonging to a specific bond.
 * \param num_bond Targeted bond index
 * \return Vector containing the indexes of the spins linked by bond num_bond.
 */
std::vector<unsigned int> Lattice1D::GetSpinsIndex(unsigned int num_bond) const{
    return std::vector<unsigned int>{num_bond,(num_bond+1)%GetBondNumber()};
}

//! Returns the number of spins in the system.
//! \return Number of spins
unsigned int Lattice1D::GetSpinNumber() const{
    return length;
}

//! Returns the number of operators in the system.
//! \return Number of bonds
unsigned int Lattice1D::GetBondNumber() const{
    return length;
}

bool Lattice1D::IsLoopStart(SSEData& ssedata, unsigned int leg) const{
    int oprtr(ssedata.GetOperator(leg/4));
    if(oprtr!=-1)
        return true;
    return false;
}

int Lattice1D::GetLoopType(SSEData& ssedata, unsigned int leg) const{
    return ssedata.GetLegSpin(leg).two_sz;
}

std::vector<double> Lattice1D::GetPathProbs(SSEData& ssedata, unsigned int leg) const{
    return std::vector<double>({.0,.0,.0,1.});
}

unsigned int Lattice1D::GetOperatorType(SSEData& ssedata, unsigned int leg) const{
    std::vector<Spin> local_vrtx(ssedata.GetBondVrtx(leg/4));
    if(local_vrtx[0].two_sz == local_vrtx[2].two_sz)
        return 0;
    return 1;
}

double Lattice1D::GetOperatorWeight(SSEData& ssedata, unsigned int leg) const{
    return 1.0;
}

bool Lattice1D::IsNegativeOperator(SSEData& ssedata, unsigned int leg) const{
    return true;
}

void SpecialUpdates(SSEData& ssedata) const{}

