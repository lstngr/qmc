#ifndef DEBUG_CONF
#define DEBUG_CONF

#include "debugconf.hpp"
#endif //DEBUG_CONF

#ifndef MCSTEP
#define MCSTEP

#include "ssedata.hpp"
#include "lattice.hpp"
#include <vector>
#include <random>

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
        protected:
        private:
            void DiagUpdate();
            void LinkVertices();
            //// Should definitely be friended to SSEData
            void OperatorLoop(std::vector<int> const& pos_first);
            void OperatorLoopSwenWang(std::vector<int> const& pos_first);
            void OperatorLoopWolff(std::vector<int> const& pos_first);
            ////
            bool AuthInsert();
            bool AuthRemove();
            
            Lattice* const lattice; //!< Lattice structure of the system. Encloses the information needed for bond index to spin indexes.
            SSEData* sse_data; //!< Spin and operator state of the system
            double beta; //!< Thermodynamic Beta
            std::string LoopMethod;
            std::mt19937 rdm_engine; //!< Random number generator engine
            std::uniform_int_distribution<unsigned int> random_bond; //!< Calls rdm_engine to generate random bonds
            std::uniform_real_distribution<double> random_prob; //!< Generates a random probability
            std::uniform_int_distribution<unsigned int> random_leg;

    }; // class MCStep

} // namespace qmc

#endif // define MCSTEP
