#include "lattice.hpp"
#include "ssedata.hpp"

namespace qmc{

    bool Spin::IsValid() const{
        if(two_norm<0)
            return false;
        else if(std::abs(two_sz) > two_norm)
            return false;
        else if(std::abs(two_sz%2) != two_norm%2)
            return false;
        return true;
    }
    void Spin::Validate(){
        if(not IsValid()){two_sz=0; two_norm=-1;}
    }
    void Spin::Raise(int amount){
        two_sz+=2*amount;
        Validate();
    }
    void Spin::Lower(int amount){
        Raise(-amount);
        Validate();
    }
    void Spin::Promote(int amount){
        two_norm+=2*amount;
        Validate();
    }
    void Spin::Demote (int amount){
        Promote(-amount);
        Validate();
    }

    bool operator==(const Spin& a, const Spin& b){
        return a.two_sz==b.two_sz && a.two_norm==b.two_norm;
    }

    std::ostream& operator<<(std::ostream& stream, const Spin& spin){
        stream << "("<< double(spin.two_norm/2) << "," << double(spin.two_sz/2.) << ")";
        return stream;
    }

    ////////////////////////////
    // LATTICE /////////////////
    ////////////////////////////

    //! Lattice class constructor
    Lattice::Lattice(){}

    ////////////////////////////
    // LATTICE1D ///////////////
    ////////////////////////////

    //! Lattice1D class constructor
    /*!
     * Instantiates a linear structure containing spins. This structure is cyclic
     * (hence, the number of spins is also the number of bonds).
     * \param length Length of the lattice, ie. the site number.
     */
    Lattice1D::Lattice1D(unsigned int length): length(length),rdm_engine(std::mt19937(5)),random_prob(std::uniform_real_distribution<double>(0.0,1.0)){}

    //! Returns spin indexes of a specific bond
    /*!
     * Returns the two indexes of each spin belonging to a specific bond.
     * \param num_bond Targeted bond index
     * \return Vector containing the indexes of the spins linked by bond num_bond.
     */
    std::vector<unsigned int> Lattice1D::GetSpinsIndex(unsigned int num_bond) const{
        return std::vector<unsigned int>{num_bond,(num_bond+1)%GetBondNumber()};
    }

    //! Returns the number of spins in the system.
    //! \return Number of spins
    unsigned int Lattice1D::GetSpinNumber() const{
        return length;
    }

    //! Returns the number of operators in the system.
    //! \return Number of bonds
    unsigned int Lattice1D::GetBondNumber() const{
        return length;
    }

    bool Lattice1D::IsLoopStart(SSEData& ssedata, unsigned int leg) const{
        int oprtr(ssedata.GetOperator(leg/4));
        if(oprtr!=-1)
            return true;
        return false;
    }

    std::vector<double> Lattice1D::FindPathProbs(SSEData& /*ssedata*/, unsigned int /*leg*/) const{
        return std::vector<double>({.0,.0,.0,1.});
    }

    int Lattice1D::FindOperatorType(SSEData& ssedata, unsigned int leg) const{
        std::vector<Spin> local_vrtx(ssedata.GetBondVrtx(leg/4));
        return FindOperatorType(local_vrtx);
    }

    int Lattice1D::FindOperatorType(std::vector<Spin> spn_state) const{
        if(spn_state.size()==2) // Append to itself if only two spins (ie diag op)
            spn_state.insert(std::end(spn_state), std::begin(spn_state), std::end(spn_state));
        if(spn_state[0]==spn_state[1]){ return -1; }
        else if(spn_state[0]==spn_state[2]){ return 0; }
        return 1;
    }

    double Lattice1D::FindOperatorWeight(SSEData& ssedata, unsigned int leg) const{
        return FindOperatorWeight(ssedata.GetOperatorType(leg/4));
    }

    double Lattice1D::FindOperatorWeight(int oprtr_type) const{
        if(oprtr_type!=-1)
            return 0.5;
        return .0;
    }

    bool Lattice1D::IsNegativeOperator(SSEData& /*ssedata*/, unsigned int /*leg*/) const{
        return false;
    }

    bool Lattice1D::IsNegativeOperator(std::vector<Spin>) const{
        return false;
    }

    void Lattice1D::SpecialUpdates(SSEData& /*ssedata*/) const{}

    bool Lattice1D::IsDiagonal(int oprtr_type) const{
        if(oprtr_type==0)
            return true;
        return false;
    }

    bool Lattice1D::IsDiagonal(SSEData& ssedata, unsigned int leg) const{
        return IsDiagonal(ssedata.GetOperatorType(leg/4));
    }

    std::vector<Spin> Lattice1D::Propagate(std::vector<Spin> spn_state, unsigned int oprtr_type) const{
        if(oprtr_type==0)
            return spn_state;
        // Else, off-diag
        if(spn_state[0].two_sz<0){
            spn_state[0].Raise(); spn_state[1].Lower();
        }else{
            spn_state[0].Lower(); spn_state[1].Raise();
        }
        return spn_state;
    }

    unsigned int Lattice1D::ChooseMove(SSEData& ssedata, unsigned int leg, Spin& loop_type){
        std::vector<double> path_probs(FindPathProbs(ssedata,leg));
        std::vector<double> comb_probs(path_probs.size());
        std::partial_sum(path_probs.begin(), path_probs.end(), comb_probs.begin());
        ssedata.SetSpinVrtx(leg,loop_type);

        // Choose Move and update
        double path_selector = random_prob(rdm_engine);
        if(path_selector < comb_probs[0]){
            // Bounce
            // Revert spin change and looptype
        }
        else if(path_selector < comb_probs[1]){
            // CS
            leg = ssedata.ContinueStraight(leg);
        }
        else if(path_selector < comb_probs[2]){
            // SWC
            leg = ssedata.SwitchAndContinue(leg);
        }
        else{
            // SWR
            leg = ssedata.SwitchAndReverse(leg);
            if(loop_type.two_sz>0)
                loop_type.Lower();
            else
                loop_type.Raise();
        }
        ssedata.SetSpinVrtx(leg,loop_type);
        return leg;
    }

    Spin Lattice1D::GetLoopType(SSEData& ssedata, unsigned int leg) const{
        Spin spn(ssedata.GetLegSpin(leg));
        if(spn.two_sz>0)
            spn.Lower();
        else
            spn.Raise();
        return spn;
    }

} // namespace qmc
