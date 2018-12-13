#include "observable.hpp"
// Other system includes
// Other project includes

using namespace qmc;

/*! Constructs a Hamiltonian operator. It mainly sets the inverse temperature,
* which is assumed constant during the lifetime of the class.
* \param beta Inverse temperature
*/
ObservableHamilton::ObservableHamilton(double beta): Observable("Energy"), beta_(beta){
    num_oprtr_ = std::vector<size_t>(); //Make sure the measurements vector is empty.
}

/*! \brief Updates the number of operators stored in
 * qmc::ObservableHamilton::num_oprtr_.
 *
 *  Pushes back the currently inserted operator number gotten from
 *  qmc::SSEData. These will be averaged on demand by
 *  qmc::ObservableHamilton::Measure.
 *  \param sim_data Pointer to the simulation data.
 */
void ObservableHamilton::UpdateObservable(SSEData* const sim_data) {
    num_oprtr_.push_back(sim_data->GetOperatorNumber());
}

/*! \brief Returns the energy of the system.
 *
 * This method averages all values stored in num_oprtr_ and returns the
 * average energy of the system.
 * \return The average value of the Hamiltonian
 */
double ObservableHamilton::Measure() const {
    double energy = - std::accumulate(num_oprtr_.begin(),
            num_oprtr_.end(), 0.0) / num_oprtr_.size() / beta_;
    return energy;
}

/* \brief Removes memorized state from the class attributes
 *
 *  It frees a bunch of memory and allows new measures (~uncorrelated) to be
 *  made.
 */
void ObservableHamilton::ClearObservable() {
    num_oprtr_.clear();
}

ObservableSign::ObservableSign()
    :Observable("Neg. Sign"), ratios(std::vector<double>()){}

void ObservableSign::UpdateObservable(SSEData* const sim_data){
    ratios.push_back(sim_data->GetNegativeWeightRatio());
}

double ObservableSign::Measure() const{
    return std::accumulate(ratios.begin(), ratios.end(), .0) / ratios.size();
}

void ObservableSign::ClearObservable() {
    ratios.clear();
}
