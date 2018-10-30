#ifndef QMC_OBSERVABLE_HPP_
#define QMC_OBSERVABLE_HPP_

#include "ssedata.hpp"

namespace qmc{

    /*! \brief Pure virtual class to instantiate OOP Observables
     * Lays the core mechanics of a physical observable.
     */
    class Observable{
        public:
            Observable(){};
            virtual ~Observable(){};
            virtual void UpdateObservable(SSEData* const sim_data)=0; /*!< Updates the attributes of the observable after each sweep. */
            virtual double Measure() const = 0; /*!< Returns the average value of the observable.*/
            virtual void ClearObservable()=0; /*!< Clears any stored variables to free space up for a new measurement.*/
    };
    
    /*! \brief Hamiltonian Observable
     * This observable evaluates the energy of the system. It stores
     * intermediate measurements locally, these require memory space increasing
     * linearly with the number of sweeps. The energy is computed as:
     * \f{equation}{
     *  \langle H \rangle = \frac{-\langle n \rangle}{\beta}
     * \f}
     * with \f$\langle n \rangle\f$ the number of operators inserted, and \f$\beta\f$
     * the inverse temperature of the system.
     */
    class ObservableHamilton : public Observable{
        public:
            ObservableHamilton(double beta);
            virtual ~ObservableHamilton(){};
            virtual void UpdateObservable(SSEData* const sim_data);
            virtual double Measure() const;
            virtual void ClearObservable();
        protected:
            std::vector<size_t> num_oprtr_; /*!<Stores the number of inserted operators for every iteration.*/
            double beta_; /*!< Inverse temperature of the system. It is set at class instantiation.*/
    };

};

#endif //QMC_OBSERVABLE_HPP_
