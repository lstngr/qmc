#ifndef MCSTEP
#define MCSTEP

#include <vector>
#include <random>

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
            std::vector<int> BondToSpins(unsigned int num_bond);
            void FlipSpins(unsigned int num_bond);
            void PrintOperators();
            void PrintSpins();
            bool AuthInsert();
            bool AuthRemove();
        protected:
        private:
            unsigned int op_str_length; //!< Limit for total operator number
            unsigned int num_cur_op; //!< Current number of operators in the system
            std::vector<int> spins; //!< Spin values
            std::vector<int> op_chain; //!< Operator string
            double beta; //!< Thermodynamic Beta

            std::mt19937 rdm_engine; //!< Random number generator engine
            std::uniform_int_distribution<unsigned int> random_bond; //!< Calls rdm_engine to generate random bonds
            std::uniform_int_distribution<unsigned int> random_spin; //!< Calls rdm_engine to generate random spins
            std::uniform_real_distribution<double> random_prob; //!< Generates a random probability

    }; // class MCStep

} // namespace qmc

#endif // define MCSTEP
