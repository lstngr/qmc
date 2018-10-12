#ifndef MCSTEP
#define MCSTEP

#include <vector>
#include <random>
#include "ssedata.hpp"

/*! \brief Class managing the QMC simulation
 *
 * This class is responsible for running Monte-Carlo sweeps on the
 * considered Heisenberg spin 1/2 system considered. It however doesn't store
 * the system in itself, which is handeld by several other classes (see
 * attributes).
 */
namespace qmc{

    class MCStep{
        public:
            MCStep();
            virtual ~MCStep();
            void DiagUpdate();
            bool AuthInsert();
            bool AuthRemove();
        protected:
        private:
            SSEData* sse_data; //!< Spin and operator state of the system
            double beta; //!< Thermodynamic Beta
            std::mt19937 rdm_engine; //!< Random number generator engine
            std::uniform_int_distribution<unsigned int> random_bond; //!< Calls rdm_engine to generate random bonds
            std::uniform_real_distribution<double> random_prob; //!< Generates a random probability

    }; // class MCStep

} // namespace qmc

#endif // define MCSTEP
