#ifndef DEBUG_CONF
#define DEBUG_CONF

#include "debug_conf.hpp"
#endif //DEBUG_CONF

#ifndef MCSTEP
#define MCSTEP

#include "ssedata.hpp"
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
            MCStep();
            virtual ~MCStep();
            void DiagUpdate();
            void LinkVertices();
            bool AuthInsert();
            bool AuthRemove();
            void SetBeta(double beta=10);
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
