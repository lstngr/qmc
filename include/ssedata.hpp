#ifndef QMC_SSEDATA_HPP_
#define QMC_SSEDATA_HPP_

#include <vector>
#include <random>
#include "lattice.hpp"
#include "debugconf.hpp"

namespace qmc{

    /*! \brief Class containing the system's spin state.
     *
     * Instances of this class are spawned by the qmc::MCStep constructor. They
     * provide storage of the simulation's state, as well as methods to display
     * the state of the system from other classes. Note that most methods
     * modifying the system are flagged private. They are meant to be called
     * from qmc::MCStep, which is thus forward-declared as friend.
     */
    class SSEData{
        public:
            SSEData(unsigned int time_slots, Lattice* const lattice, std::mt19937* const rdm_engine);
            virtual ~SSEData();
            unsigned int GetOperatorNumber(bool inserted=true) const;
            int GetSpin(unsigned int idx_spin) const; // Could combine with FlipSpins
            int GetOperator(unsigned int idx_oprtr) const;
            int GetLink(unsigned int vertex) const;
            std::vector<int> GetBondSpins(unsigned int num_bond) const;
            void PrintOperators() const;
            void PrintSpins() const;
            void PrintLinks() const;
            unsigned int GetSpinNumber() const;
            int SwitchAndReverse(int idx_leg) const;
        protected:

        private:
            void SetSpin(unsigned int idx_spin, int val_spin);
            void SetOperator(unsigned int idx_oprtr, int val_oprtr=-1);
            void FlipOperator(unsigned int idx_oprtr);
            void SetLink(unsigned int vertex, int val_vertex);
            void LinkVertex(unsigned int vertex_0, unsigned int vertex_1);
            void FlipSpins(unsigned int num_bond);
            void FlipSpin(unsigned int idx_spin);
            void ResetLinkVertex();
            void Equilibrate(double expansion_factor);

            friend class MCStep; //!< Forward declaration of qmc::MCStep. This allows moving most methods of SSEData to private, while retaining access from qmc::MCStep.
            int num_time_slots; //!< Number of time slots to place operators
            std::vector<int> spins; //!< Spin state of the system (plus or minus ones in the Heisenberg 1/2 model)
            std::vector<int> oprtr; //!< Operator string. -1 indicates the identity operator
            std::vector<int> link_vertex; //!< Linked vertex storage of the configuration.
            Lattice* const lattice; /*!< Lattice obtained from the qmc::MCStep constructor. The initialization is left to the user.*/

    };

}; // namespace qmc

#endif // #define QMC_SSEDATA_HPP_
