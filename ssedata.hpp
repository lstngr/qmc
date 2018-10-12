#ifndef SSEDATA
#define SSEDATA

#include <vector>
#include <random>

/*! \brief Class containing the system's spin state.
 *
 * A larger description is to be expected in the forthcoming developpement. Stay
 * tuned!
 */
namespace qmc{

    class SSEData{
        public:
            SSEData(unsigned int time_slots, unsigned int spin_number, std::mt19937* rdm_engine);
            virtual ~SSEData();
            unsigned int OperatorNumber(bool inserted=true);
            int GetOperator(unsigned int idx_oprtr);
            void SetOperator(unsigned int idx_oprtr, int val_oprtr=-1);
        protected:

        private:
            int num_time_slots; //!< Number of time slots to place operators
            std::vector<int> spins; //!< Spin state of the system (plus or minus ones in the Heisenberg 1/2 model)
            std::vector<int> oprtr; //!< Operator string. -1 indicates the identity operator
            std::vector<std::vector<int>> bond_vertex; //!< Linked-vertex storage of the configuration.

    };

}; // namespace qmc

#endif // #define SSEDATA
