#ifndef SSEDATA
#define SSEDATA

#include <vector>

/*! \brief Class containing the system's spin state.
 *
 * A larger description is to be expected in the forthcoming developpement. Stay
 * tuned!
 */
namespace qmc{

    class SSEData{
        public:
            SSEData();
            virtual ~SSEData();
        protected:

        private:
            int num_time_slots; //!< Number of time slots to place operators
            std::vector<int> spins; //!< Spin state of the system (plus or minus ones in the Heisenberg 1/2 model)
            std::vector<int> oprtr; //!< Operator string. -1 indicates the identity operator
            std::vector<std::vector<int>> bond_vertex; //!< Linked-vertex storage of the configuration.

    };

}; // namespace qmc

#endif // #define SSEDATA
