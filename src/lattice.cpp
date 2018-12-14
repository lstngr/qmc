#include "lattice.hpp"
#include "ssedata.hpp"

namespace qmc{

    std::uniform_real_distribution<double> Lattice::random_prob = std::uniform_real_distribution<double>(.0,1.);

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
    Spin& Spin::Raise(int amount){
        two_sz+=2*amount;
        Validate();
        return *this;
    }
    Spin& Spin::Lower(int amount){
        Raise(-amount);
        Validate();
        return *this;
    }
    Spin& Spin::Promote(int amount){
        two_norm+=2*amount;
        Validate();
        return *this;
    }
    Spin& Spin::Demote (int amount){
        Promote(-amount);
        Validate();
        return *this;
    }

    Spin& Spin::RaiseD() {
        two_sz+=2;
        if(two_norm!=0)
            Demote();
        else
            Promote();
        return *this;
    }

    Spin& Spin::LowerD() {
        two_sz-=2;
        if(two_norm!=0)
            Demote();
        else
            Promote();
        return *this;
    }

    Spin& Spin::Dz() {
        if(two_norm==0)
            Promote();
        else
            Demote();
        return *this;
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
    Lattice1D::Lattice1D(unsigned int length): length(length),rdm_engine(std::mt19937(5)){
        rdm_spin_sz = std::uniform_int_distribution<int>(0,1);
    }

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

    std::vector<double> Lattice1D::FindPathProbs(SSEData& /*ssedata*/, unsigned int /*leg*/, const Spin&) const{
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
        std::vector<double> path_probs(FindPathProbs(ssedata,leg,loop_type));
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

    Spin Lattice1D::GetLoopType(SSEData& ssedata, unsigned int leg) {
        Spin spn(ssedata.GetLegSpin(leg));
        if(spn.two_sz>0)
            spn.Lower();
        else
            spn.Raise();
        return spn;
    }

    Spin Lattice1D::RandomSpin(){
        return Spin(1,2*rdm_spin_sz(rdm_engine)-1);
    }

    ////////////////////////////
    // LATTICEDIAMOND //////////
    ////////////////////////////

    LatticeDiamond::LatticeDiamond(unsigned int length, double j2j1)
        :length(2*length), epsilon(j2j1), rdm_engine(std::mt19937(7)),
        rdm_spin_norm(std::uniform_int_distribution<unsigned int>(0,1)),
        rdm_spin_sz(std::uniform_int_distribution<int>(-1,1)){}

    std::vector<unsigned int> LatticeDiamond::GetSpinsIndex(unsigned int num_bond) const{
        return std::vector<unsigned int>({num_bond,(num_bond+1)%length});
    }

    unsigned int LatticeDiamond::GetSpinNumber() const{ return length; }
    unsigned int LatticeDiamond::GetBondNumber() const{ return length; }
    unsigned int LatticeDiamond::HowManyOperators() const { return 9; }

    bool LatticeDiamond::IsLoopStart(SSEData& ssedata, unsigned int leg) const{
        if(ssedata.GetOperatorType(leg/4)<1)
            return false; // If operator is full singlet (or empty)
        else if(ssedata.GetLegSpin(leg).two_norm < 2)
            return false; // If leg is singlet
        return true;
    }

    unsigned int LatticeDiamond::ChooseMove(SSEData& ssedata, unsigned int leg, Spin& loop_type){
        std::vector<double> path_probs(FindPathProbs(ssedata,leg,loop_type));
        std::vector<double> comb_probs(path_probs.size());
        std::partial_sum(path_probs.begin(), path_probs.end(), comb_probs.begin());
        int del_sz(loop_type.two_sz-ssedata.GetLegSpin(leg).two_sz);
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
            if(del_sz>0)
                loop_type = ssedata.GetLegSpin(leg).Raise();
            else
                loop_type = ssedata.GetLegSpin(leg).Lower();
        }
        else if(path_selector < comb_probs[2]){
            // SWC
            leg = ssedata.SwitchAndContinue(leg);
            if(del_sz>0)
                loop_type = ssedata.GetLegSpin(leg).Raise();
            else
                loop_type = ssedata.GetLegSpin(leg).Lower();
        }
        else{
            // SWR
            leg = ssedata.SwitchAndReverse(leg);
            if(del_sz>0)
                loop_type = ssedata.GetLegSpin(leg).Lower();
            else
                loop_type = ssedata.GetLegSpin(leg).Raise();
        }
        ssedata.SetSpinVrtx(leg,loop_type);
        return leg;
    }

    std::vector<double> LatticeDiamond::FindPathProbs(SSEData& ssedata, unsigned int leg, const Spin& loop_type) const{
        // XXX - Expect triplet loop type.
        std::vector<Spin> spn_state(ssedata.GetBondVrtx(leg/4));
        if(FindOperatorType(ssedata,leg)==3) // CS mandatory
            return std::vector<double>({0.,1.,0.,0.});
        // XXX - Reorder operator as AB site
        if(ssedata.GetOperatorBond(leg/4)%2){ // Location at B-site requires spin swapping
            std::iter_swap(spn_state.begin()  ,spn_state.begin()+1);
            std::iter_swap(spn_state.begin()+2,spn_state.begin()+3);
            leg = ssedata.SwitchAndReverse(leg);
        }
        unsigned int in_leg(leg%4);
        int del_sz(loop_type.two_sz-spn_state[in_leg].two_sz);
        std::vector<Spin> move_mask(spn_state);
        move_mask[in_leg] = loop_type;
        std::vector<int> out_oprtr_type(4,-1);
        // -- No Bounce --
        // -- Try CS --
        unsigned int out_leg(ssedata.ContinueStraight(in_leg));
        if(del_sz>0)
            move_mask[out_leg].Raise();
        else
            move_mask[out_leg].Lower();
        if(move_mask[out_leg].two_norm!=-1)
            out_oprtr_type[1] = FindOperatorType(move_mask);
        move_mask = spn_state;
        // -- Try SWC --
        out_leg = ssedata.SwitchAndContinue(in_leg);
        if(del_sz>0)
            move_mask[out_leg].Raise();
        else
            move_mask[out_leg].Lower();
        if(move_mask[out_leg].two_norm!=-1)
            out_oprtr_type[2] = FindOperatorType(move_mask);
        move_mask = spn_state;
        // -- Try SWR --
        out_leg = ssedata.SwitchAndReverse(in_leg);
        if(del_sz>0)
            move_mask[out_leg].Lower();
        else
            move_mask[out_leg].Raise();
        if(move_mask[out_leg].two_norm!=-1)
            out_oprtr_type[3] = FindOperatorType(move_mask);

        // Count available out operators
        unsigned int avail_traj(0);
        std::vector<double> path_probs(4,0.);
        for(unsigned int i(0);i<out_oprtr_type.size();i++){
            if(out_oprtr_type[i]==1){
                path_probs[i]=1.0;
                return path_probs;
            }else if(out_oprtr_type[i]!=-1){
                avail_traj++;
            }
        }
        for(unsigned int i(0);i<out_oprtr_type.size();i++){
            if(out_oprtr_type[i]!=-1){
                path_probs[i] = 1./avail_traj;
            }
        }
        return path_probs;
    }

    Spin LatticeDiamond::GetLoopType(SSEData& ssedata, unsigned int leg){
        Spin spn(ssedata.GetLegSpin(leg));
        if(spn.two_norm<2)
            return ssedata.GetLegSpin(leg); // FIXME - AVOID +-/-+ oprtrs
        else if(spn.two_sz>0)
            spn.Lower();
        else if(spn.two_sz<0)
            spn.Raise();
        else if(random_prob(rdm_engine)<0.5)
            spn.Raise(); // No restrictions, raise or lower.
        else
            spn.Lower();
        return spn;
    }

    int LatticeDiamond::FindOperatorType(SSEData& ssedata, unsigned int leg) const{
        std::vector<Spin> spn_state(ssedata.GetBondVrtx(leg/4));
        if(ssedata.GetOperatorBond(leg/4)%2){ // Location at B-site requires spin swapping
            std::iter_swap(spn_state.begin()  ,spn_state.begin()+1);
            std::iter_swap(spn_state.begin()+2,spn_state.begin()+3);
        }
        return FindOperatorType(spn_state);
    }

    int LatticeDiamond::FindOperatorType(std::vector<Spin> spn_state) const{
        // XXX - Assumes left leg on site A.
        if(spn_state.size()==2) // Append to itself if only two spins (ie diag op)
            spn_state.insert(std::end(spn_state), std::begin(spn_state), std::end(spn_state));
        // Check Spin Validity
        for(auto it=spn_state.begin(); it!=spn_state.end(); ++it)
            if(not (*it).IsValid())
                return -1;
        // Check Operator Validity
        if(spn_state[0].two_sz + spn_state[1].two_sz != spn_state[2].two_sz + spn_state[3].two_sz)
            return -1; // Sz conservation
        else if(std::abs(spn_state[0].two_sz - spn_state[2].two_sz) > 1 or
                std::abs(spn_state[1].two_sz - spn_state[3].two_sz) > 1)
            return -1; // Sz ladders have delta 1 or less
        else if(spn_state[1].two_norm != spn_state[3].two_norm)
            return -1; // B-Sites don't flip singlet<->triplet

        if(spn_state[0].two_norm==0 != spn_state[2].two_norm==0){ // Singlet-Triplet conversion
            if(spn_state[3].two_sz - spn_state[1].two_sz < 0)
                return 6; // DpTm
            else if(spn_state[3].two_sz - spn_state[1].two_sz > 0)
                return 7; // DmTp
            else
                return 8; // DzTz
        }
        else if(spn_state[0].two_norm==0){ // Diagonal Singlet Element
            if(spn_state[1].two_norm==0)
                return 0; // Double singlet (I16)
            else
                return 2; // I4Tz
        }
        else if(spn_state[2].two_sz - spn_state[0].two_sz > 0)
            return 4; // TpTm
        else if(spn_state[2].two_sz - spn_state[0].two_sz < 0)
            return 5; // TmTp
        else if(spn_state[0].two_sz + spn_state[1].two_sz == 0){
            return 1; // PlusMinus TzTz
        }
        // TODO - Check if more conditions needed!!
        return 3; // Remaining TzTz
    }

    double LatticeDiamond::FindOperatorWeight(SSEData& ssedata, unsigned int leg) const{
        // XXX - Negative operators are returned with positive weight. Call
        // IsNegativeOperator for sign info
        return FindOperatorWeight(FindOperatorType(ssedata,leg));
    }

    double LatticeDiamond::FindOperatorWeight(int oprtr_type) const{
        if(oprtr_type==0)
            return 1.+epsilon/2.;
        else if(oprtr_type==1)
            return epsilon;
        else if(oprtr_type==2)
            return .5*(1.+epsilon);
        else if(oprtr_type>2 and oprtr_type<7)
            return epsilon/2.;
        return 0.;
    }

    bool LatticeDiamond::IsNegativeOperator(std::vector<Spin> spn_state) const{
        // XXX - Expects AB site formatting.
        if(spn_state.size()==2) // Append to itself if only two spins (ie diag op)
            spn_state.insert(std::end(spn_state), std::begin(spn_state), std::end(spn_state));
        if(FindOperatorType(spn_state)<5)
            return false;
        else if(spn_state[0].two_sz<0 or spn_state[2].two_sz<0)
            return true; // State 5 negative operator
        else if(spn_state[1].two_sz<0 and spn_state[3].two_sz<0)
            return true; // State 6 negative operator
        return false;
    }

    bool LatticeDiamond::IsNegativeOperator(SSEData& ssedata, unsigned int leg) const{
        std::vector<Spin> spn_state(ssedata.GetBondVrtx(leg/4));
        if(ssedata.GetOperatorBond(leg/4)%2){ // Location at B-site requires spin swapping
            std::iter_swap(spn_state.begin()  ,spn_state.begin()+1);
            std::iter_swap(spn_state.begin()+2,spn_state.begin()+3);
        }
        return IsNegativeOperator(spn_state);
    }

    void LatticeDiamond::SpecialUpdates(SSEData& ssedata) const{}

    bool LatticeDiamond::IsDiagonal(SSEData& ssedata, unsigned int leg) const{
        return IsDiagonal(FindOperatorType(ssedata,leg));
    }

    bool LatticeDiamond::IsDiagonal(int oprtr_type) const{
        return oprtr_type<4;
    }

    std::vector<Spin> LatticeDiamond::Propagate(std::vector<Spin> spn_state, unsigned int oprtr_type) const{
        if(oprtr_type==4){
            spn_state[0].Raise(); spn_state[1].Lower();
        } else if(oprtr_type==5) {
            spn_state[1].Raise(); spn_state[0].Lower();
        } else if(oprtr_type==6) {
            spn_state[0].RaiseD(); spn_state[0].Lower();
        } else if(oprtr_type==7) {
            spn_state[0].LowerD(); spn_state[1].Raise();
        } else if(oprtr_type==8) {
            spn_state[0].Dz();
        }
        return spn_state;
    }

    Spin LatticeDiamond::RandomSpin() {
        unsigned int spn_norm(rdm_spin_norm(rdm_engine));
        return Spin(spn_norm,2*spn_norm*rdm_spin_sz(rdm_engine));
    }
} // namespace qmc
