//#define DEBUG_MODE
#ifndef SSEDATA
#define SSEDATA

#include <vector>
#include <random>

namespace qmc{

    /*! \brief Class containing the system's spin state.
     *
     * A larger description is to be expected in the forthcoming developpement. Stay
     * tuned!
     */
    class SSEData{
        public:
            SSEData(unsigned int time_slots, unsigned int spin_number, std::mt19937* rdm_engine);
            virtual ~SSEData();
            unsigned int GetOperatorNumber(bool inserted=true) const;
            int GetOperator(unsigned int idx_oprtr) const;
            void SetOperator(unsigned int idx_oprtr, int val_oprtr=-1);
            void LinkVertex(unsigned int vertex_0, unsigned int vertex_1);
            std::vector<int> GetBondSpins(unsigned int num_bond) const;
            void FlipSpins(unsigned int num_bond);
            void PrintOperators() const;
            void PrintSpins() const;
            void PrintLinks() const;
            unsigned int GetSpinNumber() const;
        protected:

        private:
            int num_time_slots; //!< Number of time slots to place operators
            std::vector<int> spins; //!< Spin state of the system (plus or minus ones in the Heisenberg 1/2 model)
            std::vector<int> oprtr; //!< Operator string. -1 indicates the identity operator
            std::vector<int> link_vertex; //!< Linked vertex storage of the configuration.

    };

}; // namespace qmc

#endif // #define SSEDATA
