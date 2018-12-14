#ifndef QMC_OBSERVABLE_HPP_
#define QMC_OBSERVABLE_HPP_

#include "ssedata.hpp"

namespace qmc{

    /*! \brief Pure virtual class to instantiate OOP Observables
     *
     * Lays the core mechanics of a physical observable.
     */
    class Observable{
        public:
            //! Base class constructor needs to set its name (as private
            //! attribute)
            //! \param name Name of the observable to print it later.
            Observable(std::string name=""):sign(std::vector<short>({})),obsrv_name_(name){};
            virtual ~Observable(){};
            virtual void UpdateObservable(SSEData* const sim_data, short sign)=0; /*!< Updates the attributes of the observable after each sweep.
            \param sim_data Pointer to a simulation container. Allows access to simulation data.*/
            virtual double Measure() const = 0; /*!< Returns the average value of the observable.
            \return The average value of the observable*/

            virtual void ClearObservable(); /*!< Clears any stored variables to free space up for a new measurement.*/
            //! Returns the name of an observable
            //! \return Name of the observable
            std::string GetName() const {return obsrv_name_;}
            void AddSign(short is_neg=false);
            double MeanSign() const;
        private:
            std::vector<short> sign;
            std::string obsrv_name_; //!< Name of the observable. Used when printing
    };

    /*! \brief Hamiltonian Observable
     *
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
            virtual void UpdateObservable(SSEData* const sim_data,short sign);
            virtual double Measure() const;
            virtual void ClearObservable();
        protected:
            std::vector<size_t> num_oprtr_; /*!<Stores the number of inserted operators for every iteration.*/
            double beta_; /*!< Inverse temperature of the system. It is set at class instantiation.*/
    };

    class ObservableSign : public Observable{
        public:
            ObservableSign();
            virtual ~ObservableSign(){};
            void UpdateObservable(SSEData* const sim_data, short sign);
            double Measure() const;
            void ClearObservable();
        protected:
        private:
            std::vector<double> ratios;
    };

};

#endif //QMC_OBSERVABLE_HPP_
