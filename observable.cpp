#include "observable.hpp"
// Other system includes
// Other project includes

using namespace qmc;

ObservableHamilton::ObservableHamilton(double beta): beta_(beta){}

void ObservableHamilton::UpdateObservable(SSEData* const sim_data) {
    num_oprtr_.push_back(sim_data->GetOperatorNumber());
}

double ObservableHamilton::Measure() const {
    double energy = - std::accumulate(num_oprtr_.begin(), 
            num_oprtr_.end(), 0.0) / num_oprtr_.size() / beta_;
    return energy;
}

