#ifndef QMC_SSEDATA_HPP_
#define QMC_SSEDATA_HPP_

#include <vector>
#include <random>
#include "debugconf.hpp"

namespace qmc{
    
    class Lattice;
    struct Spin;

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
            Spin GetSpin(unsigned int idx_spin) const;
            void SetSpin(unsigned int idx_spin, int two_s, int two_sz);
            void SetSpin(unsigned int idx_spin, Spin spn);
            void SetSpinVrtx(unsigned int idx_leg, Spin spn);
            Spin GetLegSpin(unsigned int leg_idx) const;
            std::vector<Spin> GetBondVrtx(unsigned int oprtr_idx) const;
            int GetOperator(unsigned int idx_oprtr) const;
            int GetOperatorType(unsigned int idx_oprtr) const;
            int GetOperatorBond(unsigned int idx_oprtr) const;
            int GetLink(unsigned int vertex) const;
            std::vector<Spin> GetBondSpins(unsigned int num_bond) const;
            void PrintOperators() const;
            void PrintSpins() const;
            void PrintLinks() const;
            unsigned int SwitchAndReverse (unsigned int idx_leg) const;
            unsigned int ContinueStraight (unsigned int idx_leg) const;
            unsigned int SwitchAndContinue(unsigned int idx_leg) const;

            void SetOperator(unsigned int idx_oprtr, int oprtr_type=-1, int oprtr_bond=0, bool update_bond_vrtx=false);
            void FlipOperator(unsigned int idx_oprtr); // Deprecated
            void SetLink(unsigned int vertex, int val_vertex);
            void LinkVertex(unsigned int vertex_0, unsigned int vertex_1);
            void UpdateSpins(unsigned int idx_oprtr);
            void RaiseSpin(unsigned int idx_leg, int amount=1);
            void LowerSpin(unsigned int idx_leg, unsigned int amount=1);
            void ResetLinkVertex();
            int  GetVFirst(unsigned int idx_spin) const;
            void SetVFirst(unsigned int idx_spin, unsigned int idx_leg);
            int  GetVLast (unsigned int idx_spin) const;
            void SetVLast(unsigned int idx_spin, unsigned int idx_leg);
            void UpdateSpinState();
            void Equilibrate(double expansion_factor);
            bool CheckPositiveWeight();

        private:
            int num_time_slots; //!< Number of time slots to place operators
            int num_oprtr_type_; //!< Number of allowed operator types
            std::vector<Spin> spins; //!< Spin state of the system (plus or minus ones in the Heisenberg 1/2 model)
            std::vector<int> oprtr; //!< Operator string. -1 indicates the identity operator
            std::vector<int> link_vertex; //!< Linked vertex storage of the configuration.
            std::vector<std::vector<Spin>> bond_vertex_; /*!< Bond vertex of the system.
            We call bond vertex an entity storing the four spins interacting at the legs of every operator.*/
            std::vector<int> v_first_; //!< Container recording the first leg at which a spin interacts. If no interaction happens, the value -1 is assigned.
            std::vector<int> v_last_; //!< Container recording the last leg at which a spin interacts. If no interaction happens, the value -1 is assigned.
            Lattice* const lattice; /*!< Lattice obtained from the qmc::MCStep constructor. The initialization is left to the user.*/

    };

}; // namespace qmc

#endif // #define QMC_SSEDATA_HPP_
