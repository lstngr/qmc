#ifndef QMC_MCSTEP_HPP_
#define QMC_MCSTEP_HPP_

#include <vector>
#include <random>
#include <ostream>
#include "debugconf.hpp"
#include "lattice.hpp"
#include "ssedata.hpp"
#include "observable.hpp"
#include "omanager.hpp"

namespace qmc{

    /*! \brief Class managing the QMC simulation
     *
     * This class is responsible for running Monte-Carlo sweeps on the
     * considered Heisenberg spin 1/2 system considered. It however doesn't store
     * the system in itself, which is handeld by several other classes (see
     * attributes).
     */
    class MCStep{
        public:
            MCStep(Lattice* const lattice);
            virtual ~MCStep();
            void Sweep(unsigned int num_sweep=1);
            void SetBeta(double beta=10);
            void AddObservable(Observable* obsrv);
            double GetBeta() const;
            const std::vector<Observable*>& GetObservables() const;
        protected:
        private:
            void DiagUpdate();
            void LinkVertices();
            //// Should definitely be friended to SSEData
            void OperatorLoop(std::vector<int> const& pos_first, bool all_links_empty=false);
            void OperatorLoopSwenWang(std::vector<int> const& pos_first);
            void OperatorLoopWolff(std::vector<int> const& pos_first, bool all_links_empty=false);
            ////
            bool AuthInsert();
            bool AuthRemove();
            std::vector<unsigned int> LegToSpins(unsigned int idx_leg) const;
            void UpdateObservables() const;

            Lattice* const lattice; //!< Lattice structure of the system. Encloses the information needed for bond index to spin indexes.
            SSEData* sse_data; //!< Spin and operator state of the system
            OManager* output_manager_;
            double beta; //!< Thermodynamic Beta
            double expansion_factor; //!< Controls the operator length
            std::string LoopMethod; /*!< Stores the OperatorLoop method to use.
                                      In its current state, the program supports:
                                      - "swenwang": Swendsen-Wang Algorithm
                                      - "wolff": Wolff Algorithm (default)*/
            std::mt19937 rdm_engine; //!< Random number generator engine
            std::uniform_int_distribution<unsigned int> random_bond; //!< Calls rdm_engine to generate random bonds
            std::uniform_real_distribution<double> random_prob; //!< Generates a random probability
            std::uniform_int_distribution<unsigned int> random_leg; //!< Generates a random bond index (for qmc::MCStep::OperatorLoopWolff)
            std::vector<Observable*> observables_; /*!< Contains pointers to user defined observables. They will be updated for every sweep, so that the user can request their values through the appropriate method afterwards.*/

    }; // class MCStep

} // namespace qmc

#endif // define QMC_MCSTEP_HPP_
