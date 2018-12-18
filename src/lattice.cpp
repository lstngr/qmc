#include "lattice.hpp"
#include "ssedata.hpp"

namespace qmc{

    std::uniform_real_distribution<double> Lattice::random_prob = std::uniform_real_distribution<double>(.0,1.);
    std::mt19937 Lattice::rdm_engine = std::mt19937(7);

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

    unsigned int Lattice::ChooseMove(SSEData& ssedata, unsigned int leg, unsigned int leg_start, Spin& loop_type){
        int del_sz(loop_type.two_sz-ssedata.GetLegSpin(leg).two_sz);
        if(leg==leg_start){
            if(del_sz==0)
                return leg_start;
            else
                throw std::logic_error(std::string("Loop is closed with non-zero delta sz."));
        }
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
    ////////////////////////////
    // LATTICE1D ///////////////
    ////////////////////////////

    //! Lattice1D class constructor
    /*!
     * Instantiates a linear structure containing spins. This structure is cyclic
     * (hence, the number of spins is also the number of bonds).
     * \param length Length of the lattice, ie. the site number.
     */
    Lattice1D::Lattice1D(unsigned int length): length(length){
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

    int Lattice1D::FindBondType(SSEData& ssedata, unsigned int bond) const{
        std::vector<Spin> spn_state(ssedata.GetBondSpins(bond));
        return FindOperatorType(spn_state);
    }

    std::vector<Spin> Lattice1D::PropagateBond(SSEData& ssedata, unsigned int oprtr) const{
        return Propagate(
                ssedata.GetBondSpins(ssedata.GetOperatorBond(oprtr)),
                ssedata.GetOperatorType(oprtr)
                );
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

    void Lattice1D::SpecialUpdates(SSEData& /*ssedata*/) {}

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
        :length(2*length), epsilon(j2j1),
        rdm_spin_norm(std::uniform_int_distribution<unsigned int>(0,3)),
        rdm_spin_sz(std::uniform_int_distribution<int>(0,2)){}

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

    std::vector<double> LatticeDiamond::FindPathProbs(SSEData& ssedata, unsigned int leg, const Spin& loop_type) const{
        if(loop_type.two_norm<2)
            throw std::logic_error(std::string("Expected triplet loop_type."));
        std::vector<Spin> spn_state(ssedata.GetBondVrtx(leg/4));
        if(FindOperatorType(ssedata,leg)==2) // CS mandatory
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
        move_mask[in_leg] = loop_type;
        // -- Try SWC --
        out_leg = ssedata.SwitchAndContinue(in_leg);
        if(del_sz>0)
            move_mask[out_leg].Raise();
        else
            move_mask[out_leg].Lower();
        if(move_mask[out_leg].two_norm!=-1)
            out_oprtr_type[2] = FindOperatorType(move_mask);
        move_mask = spn_state;
        move_mask[in_leg] = loop_type;
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
            throw std::logic_error(std::string("Shouldn't ask loop_type on singlet legs"));
        else if(spn.two_sz>0)
            spn.Lower();
        else if(spn.two_sz<0)
            spn.Raise();
        else if(random_prob(rdm_engine)<0.5)
            // Current spin is null triplet
            spn.Raise();
        else
            spn.Lower();
        return spn;
    }

    std::vector<Spin> LatticeDiamond::PropagateBond(SSEData& ssedata, unsigned int oprtr) const{
        unsigned int bond(ssedata.GetOperatorBond(oprtr));
        std::vector<Spin> spn_state(ssedata.GetBondSpins(bond));
        if(bond%2){ // Location at B-site requires spin swapping
            std::iter_swap(spn_state.begin()  ,spn_state.begin()+1);
        }
        spn_state = Propagate(spn_state,ssedata.GetOperatorType(oprtr));
        if(bond%2){ // Location at B-site requires spin swapping
            std::iter_swap(spn_state.begin()  ,spn_state.begin()+1);
        }
        return spn_state;
    }

    int LatticeDiamond::FindBondType(SSEData& ssedata, unsigned int bond) const{
        std::vector<Spin> spn_state(ssedata.GetBondSpins(bond));
        if(bond%2) // Location at B-site requires spin swapping
            std::iter_swap(spn_state.begin(),spn_state.begin()+1);
        return FindOperatorType(spn_state);
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
        else if(std::abs(spn_state[0].two_sz - spn_state[2].two_sz) > 2 or
                std::abs(spn_state[1].two_sz - spn_state[3].two_sz) > 2)
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
        else if(spn_state[0].two_norm==0 or spn_state[1].two_norm==0){ // Diagonal Singlet Element
            if(spn_state[0].two_norm==0 and spn_state[1].two_norm==0)
                return 0; // Double singlet (I16)
            else
                return 2; // I4Tz or TzI4
        }
        else if(spn_state[2].two_sz - spn_state[0].two_sz > 0)
            return 4; // TpTm
        else if(spn_state[2].two_sz - spn_state[0].two_sz < 0)
            return 5; // TmTp
        else if(spn_state[0].two_sz + spn_state[1].two_sz == 0
                and spn_state[0].two_sz != 0){
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
        else if(oprtr_type>2 and oprtr_type<9)
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

    void LatticeDiamond::SpecialUpdates(SSEData& ssedata) {
        ASiteUpdate(ssedata);
        BSiteUpdate(ssedata);
    }

    void LatticeDiamond::ASiteUpdate(SSEData& ssedata) {
        for(unsigned int idx_spin(0); idx_spin<this->GetSpinNumber();
                idx_spin+=2){
            int init_leg(ssedata.GetVFirst(idx_spin));
            int clst_leg(init_leg); // Cluster taging leg
            double prob(1.);
            if(init_leg!=-1){ // Otherwise updated in Sweep
                std::vector<int> old_oprtr_type({});
                std::vector<int> new_oprtr_type({});
                std::vector<unsigned int> clst_legs({});
                // Fill type difference vectors
                do{
                    int some_type(ssedata.GetOperatorType(clst_leg/4));
                    if(some_type==1) some_type=-1;
                    else if((some_type==2 or some_type==3) and ssedata.GetLegSpin(clst_leg).two_sz!=0) some_type=-1;
                    old_oprtr_type.push_back(some_type);
                    switch(some_type){
                        case 0: some_type=2; break;
                        case 3: some_type=2; break;
                        case 4: some_type=6; break;
                        case 5: some_type=7; break;
                        case 6: some_type=4; break;
                        case 7: some_type=5; break;
                        case 2:
                                std::vector<Spin> spn_state(ssedata.GetBondVrtx(clst_leg/4));
                                if(spn_state[clst_leg%4].two_norm==0){ //Already know spin is in cluster
                                    spn_state[clst_leg%4].Promote();
                                    spn_state[ssedata.ContinueStraight(clst_leg)%4].Promote();
                                } else {
                                    spn_state[clst_leg%4].Demote();
                                    spn_state[ssedata.ContinueStraight(clst_leg)%4].Demote();
                                }
                                if(clst_leg%2){ // Operator acts on BA site
                                    std::iter_swap(spn_state.begin()  ,spn_state.begin()+1);
                                    std::iter_swap(spn_state.begin()+2,spn_state.begin()+3);
                                }
                                some_type = this->FindOperatorType(spn_state);
                                break;
                    }
                    new_oprtr_type.push_back(some_type);
                    clst_legs.push_back(clst_leg);
                    clst_leg = ssedata.GetLink(ssedata.ContinueStraight(clst_leg));
                }while(clst_leg!=init_leg);

                // Check something can actually be flipped
                bool may_flip(false);
                for(unsigned int idx_oprtr(0); idx_oprtr<old_oprtr_type.size(); idx_oprtr++){
                    if(old_oprtr_type[idx_oprtr]!=-1){may_flip=true; break;}
                }

                // Compute cluster flip probabilities.
                bool clst_computed(false);
                int first_clst(-1);
                int idx_oprtr(0);
                unsigned int num_oprtr(old_oprtr_type.size());
                do{
                    unsigned int left_idx ((idx_oprtr+num_oprtr-1)%num_oprtr);
                    unsigned int right_idx((idx_oprtr+num_oprtr+1)%num_oprtr);
                    if(not may_flip) break;
                    else if(new_oprtr_type[idx_oprtr]==-1) idx_oprtr = right_idx;
                    else if(new_oprtr_type[left_idx]!=-1){
                        idx_oprtr = left_idx;
                        if(idx_oprtr==0){
                            // We browsed the whole set without finding edges,
                            // compute the whole flip probability.
                            for(unsigned int idx_flip(0); idx_flip<num_oprtr; idx_flip++){
                                prob *= this->FindOperatorWeight(new_oprtr_type[idx_flip]) / this->FindOperatorWeight(old_oprtr_type[idx_flip]);
                            }
                            if(random_prob(rdm_engine)<prob){
                                for(unsigned int idx_flip(0); idx_flip<num_oprtr; idx_flip++){
                                    // Iterate through cluster data and flip spins
                                    ssedata.SetOperator(clst_legs[idx_flip]/4, new_oprtr_type[idx_flip],
                                            ssedata.GetOperatorBond(clst_legs[idx_flip]/4),false);
                                    Spin leg_spin(ssedata.GetLegSpin(clst_legs[idx_flip]));
                                    if(leg_spin.two_sz==0){
                                        if(leg_spin.two_norm==0)
                                            ssedata.SetSpinVrtx(clst_legs[idx_flip],leg_spin.Promote());
                                        else
                                            ssedata.SetSpinVrtx(clst_legs[idx_flip],leg_spin.Demote());
                                    }
                                    leg_spin=ssedata.GetLegSpin(ssedata.ContinueStraight(clst_legs[idx_flip]));
                                    if(leg_spin.two_sz==0){
                                        if(leg_spin.two_norm==0)
                                            ssedata.SetSpinVrtx(ssedata.ContinueStraight(clst_legs[idx_flip]),leg_spin.Promote());
                                        else
                                            ssedata.SetSpinVrtx(ssedata.ContinueStraight(clst_legs[idx_flip]),leg_spin.Demote());
                                    }
                                }
                            }
                            clst_computed = true;
                        }
                    }
                    else if(idx_oprtr==first_clst) clst_computed=true; // We are on the edge of the first cluster we found. Exit loop
                    else{ // We sit on (bottom) cluster edge
                        if(first_clst==-1) //Record first cluster encounter
                            first_clst = idx_oprtr;
                        unsigned int clst_start(idx_oprtr);
                        do{ // Compute flip probability
                            prob *= this->FindOperatorWeight(new_oprtr_type[idx_oprtr]) / this->FindOperatorWeight(old_oprtr_type[idx_oprtr]);
                            idx_oprtr = (idx_oprtr+num_oprtr+1)%num_oprtr; // Move right
                        }while(new_oprtr_type[idx_oprtr]==-1); // We reach end of cluster
                        if(random_prob(rdm_engine)<prob){ // Decide flipping the cluster
                            for(int idx_flip(clst_start); idx_flip<idx_oprtr; idx_flip++){
                                // Iterate through cluster data and flip spins
                                ssedata.SetOperator(clst_legs[idx_flip]/4, new_oprtr_type[idx_flip],
                                        ssedata.GetOperatorBond(clst_legs[idx_flip]/4),false);
                                Spin leg_spin(ssedata.GetLegSpin(clst_legs[idx_flip]));
                                if(leg_spin.two_norm==0)
                                    ssedata.SetSpinVrtx(clst_legs[idx_flip],leg_spin.Promote());
                                else
                                    ssedata.SetSpinVrtx(clst_legs[idx_flip],leg_spin.Demote());
                                leg_spin=ssedata.GetLegSpin(ssedata.ContinueStraight(clst_legs[idx_flip]));
                                if(leg_spin.two_norm==0)
                                    ssedata.SetSpinVrtx(ssedata.ContinueStraight(clst_legs[idx_flip]),leg_spin.Promote());
                                else
                                    ssedata.SetSpinVrtx(ssedata.ContinueStraight(clst_legs[idx_flip]),leg_spin.Demote());
                            }
                        }
                        prob = 1.;
                    }
                }while(not clst_computed);
            }
        }
    }

    void LatticeDiamond::BSiteUpdate(SSEData& ssedata) {
        for(unsigned int idx_spin(1);
                idx_spin<this->GetSpinNumber();idx_spin+=2){
            int init_leg(ssedata.GetVFirst(idx_spin));
            double prob(1.);
            if(init_leg!=-1) {// Spin flipping is taken care of by DiagUpdate
                // Initialize "loop type"
                Spin curSpin(ssedata.GetLegSpin(init_leg));
                Spin loop_type(Spin(-1,0));
                if(curSpin.two_norm==0)
                    loop_type=Spin(2,0);
                else
                    loop_type=Spin(0,0);

                // Compute Flip Probability
                int leg_idx(init_leg);
                bool may_flip(true);
                do{
                    int oprtr_type(ssedata.GetOperatorType(leg_idx/4));
                    if(not this->IsDiagonal(oprtr_type)){
                        may_flip = false;
                        break; // Spin is changing, abort
                    } else if(oprtr_type!=2)
                        prob *= this->FindOperatorWeight(2)/this->FindOperatorWeight(oprtr_type);
                    else{
                        std::vector<Spin> spn_state(ssedata.GetBondVrtx(leg_idx/4));
                        int in_leg(leg_idx%4);
                        spn_state[in_leg] = loop_type;
                        spn_state[ssedata.ContinueStraight(in_leg)] = loop_type;
                        prob *= this->FindOperatorWeight(this->FindOperatorType(spn_state)) / this->FindOperatorWeight(2);
                    }
                    leg_idx = ssedata.GetLink(ssedata.ContinueStraight(leg_idx));
                }while(leg_idx!=init_leg);

                // Try flip
                double flip_prob(random_prob(rdm_engine));
                if(may_flip and flip_prob<prob){
                    do{
                        ssedata.SetSpinVrtx(leg_idx,loop_type);
                        ssedata.SetSpinVrtx(ssedata.ContinueStraight(leg_idx),loop_type);
                        ssedata.SetOperator(leg_idx/4,this->FindOperatorType(ssedata,leg_idx),
                                ssedata.GetOperatorBond(leg_idx/4),false);
                        leg_idx = ssedata.GetLink(ssedata.ContinueStraight(leg_idx));
                    }while(leg_idx!=init_leg);
                }
            }
        }
    }

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
        unsigned int spn_norm(2*(rdm_spin_norm(rdm_engine)!=0));
        return Spin(spn_norm,spn_norm*(rdm_spin_sz(rdm_engine)-1));
    }
} // namespace qmc
