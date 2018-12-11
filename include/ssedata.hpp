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
            int SwitchAndReverse (int idx_leg) const;
            int ContinueStraight (int idx_leg) const;
            int SwitchAndContinue(int idx_leg) const;
        
            void SetSpin(unsigned int idx_spin, int val_spin); // To move to s,sz params
            void SetSpin(Spin spin_var);
            void SetSz(int sz);
            void SetSNorm(int s_norm);
            void SetOperator(unsigned int idx_oprtr, int val_oprtr=-1);
            void FlipOperator(unsigned int idx_oprtr); // Deprecated
            void SetLink(unsigned int vertex, int val_vertex);
            void LinkVertex(unsigned int vertex_0, unsigned int vertex_1);
            void FlipSpins(unsigned int num_bond); // Deprecated
            void FlipSpin(unsigned int idx_spin); // Deprecated
            void RaiseSpin(unsigned int idx_spin, unsigned int amount=1);
            void LowerSpin(unsigned int idx_spin, unsigned int amount=1);
            void ResetLinkVertex();
            void Equilibrate(double expansion_factor);

        private:
            int num_time_slots; //!< Number of time slots to place operators
            int num_oprtr_type_;
            std::vector<Spin> spins; //!< Spin state of the system (plus or minus ones in the Heisenberg 1/2 model)
            std::vector<int> oprtr; //!< Operator string. -1 indicates the identity operator
            std::vector<int> link_vertex; //!< Linked vertex storage of the configuration.
            std::vector<std::vector<Spin>> bond_vertex_;
            std::vector<unsigned int> v_first_;
            std::vector<unsigned int> v_last_;
            Lattice* const lattice; /*!< Lattice obtained from the qmc::MCStep constructor. The initialization is left to the user.*/

    };

}; // namespace qmc

#endif // #define QMC_SSEDATA_HPP_
