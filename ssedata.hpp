#ifndef DEBUG_CONF
#define DEBUG_CONF
#include "debugconf.hpp"
#endif //DEBUG_CONF

#ifndef SSEDATA
#define SSEDATA

#include <vector>
#include <random>
#include "lattice.hpp"

namespace qmc{

    /*! \brief Class containing the system's spin state.
     *
     * A larger description is to be expected in the forthcoming developpement. Stay
     * tuned!
     */
    class SSEData{
        public:
            SSEData(unsigned int time_slots, Lattice* const lattice, std::mt19937* const rdm_engine);
            virtual ~SSEData();
            unsigned int GetOperatorNumber(bool inserted=true) const;
            //// Could combine with FlipSpins
            int GetSpin(unsigned int idx_spin) const;
            void SetSpin(unsigned int idx_spin, int val_spin);
            ////
            int GetOperator(unsigned int idx_oprtr) const;
            //// Could maybe friend?
            void SetOperator(unsigned int idx_oprtr, int val_oprtr=-1);
            void FlipOperator(unsigned int idx_oprtr);
            ////
            int GetLink(unsigned int vertex) const;
            void SetLink(unsigned int vertex, int val_vertex);
            void LinkVertex(unsigned int vertex_0, unsigned int vertex_1);
            std::vector<int> GetBondSpins(unsigned int num_bond) const;
            void FlipSpins(unsigned int num_bond);
            void FlipSpin(unsigned int idx_spin);
            void PrintOperators() const;
            void PrintSpins() const;
            void PrintLinks() const;
            unsigned int GetSpinNumber() const;
            void ResetLinkVertex();
            int SwitchAndReverse(int idx_leg) const;
            void Equilibrate(double expansion_factor);
        protected:

        private:
            int num_time_slots; //!< Number of time slots to place operators
            std::vector<int> spins; //!< Spin state of the system (plus or minus ones in the Heisenberg 1/2 model)
            std::vector<int> oprtr; //!< Operator string. -1 indicates the identity operator
            std::vector<int> link_vertex; //!< Linked vertex storage of the configuration.
            Lattice* const lattice; /*!< Lattice obtained from the qmc::MCStep constructor. The initialization is left to the user.*/

    };

}; // namespace qmc

#endif // #define SSEDATA
