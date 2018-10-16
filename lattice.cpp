#include "lattice.hpp"

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

// Lattice1D class constructor
Lattice1D::Lattice1D(){

}

//! Returns spin indexes of a specific bond
/*!
 * Returns the two indexes of each spin belonging to a specific bond.
 * \param num_bond Targeted bond index
 */
std::vector<unsigned int> Lattice1D::GetSpinsIndex(unsigned int num_bond) const{
    return std::vector<unsigned int>{num_bond,num_bond+1};
}
