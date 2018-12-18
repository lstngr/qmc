#ifndef QMC_LATTICE_HPP_
#define QMC_LATTICE_HPP_

#include <vector>
#include <ostream>
#include <random>

namespace qmc{

    class SSEData;

    /*! Memory representation of a spin. It allows a more general handling than
     * a simple vector, as well as the ability to assign (C++) operators to the
     * class to simplify th program's syntax. Furthermore, it provides a way to
     * quickly flag invalid spins.*/
    struct Spin{
        //! Builds a qmc::Spin object. Note one must pass the TWO-norm and
        //! TWO-z-component of the spin to this function. This convention is
        //! simply due to the fact that using integers reduces memory use (and
        //! numerical error) over doubles.
        Spin(int tsn, int tsz):two_norm(tsn),two_sz(tsz){}
        Spin& Raise(int amount=1);
        Spin& Lower(int amount=1);
        Spin& Promote(int amount=1);
        Spin& Demote (int amount=1);
        Spin& RaiseD ();
        Spin& LowerD ();
        Spin& Dz();
        void Validate();
        bool IsValid() const;
        int two_norm; //!< Norm of the spin
        int two_sz; //!< Z-Component of the spin
    };
    bool operator==(const Spin& a, const Spin& b);
    std::ostream& operator<<(std::ostream& stream, const Spin& spin);

    //! Base class managing Lattice structures.
    /*!
     * This class is meant to enable OOP programmation later in the code. It
     * will for example be called in the initalization, and can provide:
     * - Consitent Bond<->Spin indexes conversion
     * - Return consistent dimensions when allocating space for SSEData objects.
     *
     * Some methods are consequently defined as pure virtual ones.
     * It is also to mention that, due to implementation restrictions in
     * qmc::MCStep, all the methods of this class will be flagged const. This
     * reflects the fact that the Lattice structure shouldn't change during a
     * running Monte-Carlo simulation.
     */
    class Lattice{
        public:
            Lattice();
            virtual ~Lattice(){};

            //! Returns spin indexes connected by num_bond
            /*!
             * For a given bond, a std::vector containing the spin indexes
             * belonging to the specified bond is returned. This function is
             * adapted polymorphically according to the lattice type.
             * \param num_bond Index of the bond
             * \return Vector containing spin indexes connected by num_bond
             */
            virtual std::vector<unsigned int> GetSpinsIndex(unsigned int num_bond) const = 0;

            //! Returns the total number of spins
            /*!
             * Returns the number of spins (or sites) enclosed in the Lattice
             * structure. This function is adapted polymorphically according to
             * the lattice type.
             * \return Number of spins
             */
            virtual unsigned int GetSpinNumber() const = 0;

            //! Returns the total number of bonds
            /*!
             * Returns the number of bonds (linking sites) enclosed in the
             * Lattice structure. This function is adapted polymorphically
             * according to the lattice type.
             * \return Number of bonds between spins.
             */
            virtual unsigned int GetBondNumber() const = 0;

            //! Returns the number of operator types available. Useful to
            //! initialize some qmc::SSEData attributes.
            virtual unsigned int HowManyOperators() const = 0;

            //! Checks if the given leg can be used as starting position of a
            //! loop
            virtual bool IsLoopStart(SSEData& ssedata, unsigned int leg) const = 0;

            //! Returns a qmc::Spin to be used by the caller loop algorithm
            //! function.
            virtual Spin GetLoopType(SSEData& ssedata, unsigned int leg) = 0;

            /*! \brief Picks an out-leg given the entry-leg and loop type.
             * Returns the exit leg when a loop (carrying spin loop_type)
             * encounters the operator located at a given leg. The loop starting
             * leg is also provided in order to account for special loop-closing
             * criterion.
             */
            virtual unsigned int ChooseMove(SSEData& ssedata, unsigned int leg, unsigned int leg_start, Spin& loop_type);

            //! Returns the probabilities of moves when encountering an operator
            //! in the following order: (bounce, straight-and-continue,
            //! switch-and-continue, switch-and-reverse).
            virtual std::vector<double> FindPathProbs(SSEData& ssedata, unsigned int leg, const Spin& loop_type) const = 0;

            //! Returns the operator type located at a given leg. This function
            //! analyses the spins attached to the specified operator to return
            //! the type. It DOESN'T fetch the type from qmc::SSEData!
            virtual int FindOperatorType(SSEData& ssedata, unsigned int leg) const = 0;

            //! Returns the type of diagonal operator that could be inserted
            //! between spins connected by the given bond. This function requires
            //! an up-to-date spin state to be stored in qmc::SSEData and meant
            //! to be called from qmc::MCStep::DiagUpdate.
            virtual int FindBondType(SSEData& ssedata, unsigned int bond) const = 0;

            //! Given an operator position, fetches the spins attached to the
            //! former operator, compute their counterparts after the operator
            //! has acted upon them, and returns a size 2 vector containing
            //! these.
            virtual std::vector<Spin> PropagateBond(SSEData& ssedata, unsigned int oprtr) const = 0;

            /*!Finds the type of the requested operator through
             * qmc::Lattice::FindOperatorType and calls
             * qmc::Lattice::FindOperatorWeight(int) to get its weight.*/
            virtual double FindOperatorWeight(SSEData& ssedata, unsigned int leg) const = 0;

            /*! Given an operator type, return its weight (simple switch-case).
             * All the weights returned by this function are positive.*/
            virtual double FindOperatorWeight(int oprtr_type) const = 0;

            /*! Given a leg, flags an operator with negative weight.*/
            virtual bool IsNegativeOperator(SSEData& ssedata, unsigned int leg) const = 0;

            /*! If special update mechanisms are required, they should be
             * inserted here. This function is called after the
             * qmc::MCStep::OperatorLoop and before updating the spin state or
             * measuring of the observables.*/
            virtual void SpecialUpdates(SSEData& ssedata) = 0;

            /*! Checks if a given operator is diagonal.*/
            virtual bool IsDiagonal(SSEData& ssedata, unsigned int leg) const = 0;

            //! Checks if a given operator type is diagonal.
            virtual bool IsDiagonal(int oprtr_type) const = 0;

            //! Returns a random spin allowed within the lattice considered.
            virtual Spin RandomSpin() = 0;
        protected:
            /*! Flags operators associated with negative weights. This function
             * is protected as it is passed a std::vector<qmc::Spin> directly.*/
            virtual bool IsNegativeOperator(std::vector<Spin> spn_state) const = 0;

            /*! Finds the weight associated to a certain operator by analyizing
             * its bond vertex. This function is protected as it is passed a
             * std::vector<qmc::Spin> directly.*/
            virtual int FindOperatorType(std::vector<Spin> spn_state) const = 0;

            /*! Computes the propagation of a pair of spins through a given
             * operator type. This function is protected as it is passed a
             * std::vector<qmc::Spin> directly.*/
            virtual std::vector<Spin> Propagate(std::vector<Spin> spn_state, unsigned int oprtr_type) const = 0;
            static std::mt19937 rdm_engine; /*!< Random number generator engine. It is shared by all qmc::Lattice instances.*/
            static std::uniform_real_distribution<double>random_prob; /*!< Uniform random generator of real numbers between 0 and 1. It is used for choosing paths and providing random spins. This object is shared by all instances of base class qmc::Lattice.*/
        private:
    };

    //! Implementation of a linear one dimensionnal lattice structure.
    /*!
     * Name is pretty self-explanatory
     */
    class Lattice1D : public Lattice{
        public:
            Lattice1D(unsigned int length);
            virtual ~Lattice1D(){};
            std::vector<unsigned int> GetSpinsIndex(unsigned int num_bond) const override;
            unsigned int GetSpinNumber() const override;
            unsigned int GetBondNumber() const override;
            unsigned int HowManyOperators() const override{return 2;}
            bool IsLoopStart(SSEData& ssedata, unsigned int leg) const override;
            Spin GetLoopType(SSEData& ssedata, unsigned int leg) override;
            std::vector<double> FindPathProbs(SSEData& ssedata, unsigned int leg, const Spin& loop_type) const override;
            int FindOperatorType(SSEData& ssedata, unsigned int leg) const override;
            int FindBondType(SSEData& ssedata, unsigned int bond) const override;
            std::vector<Spin> PropagateBond(SSEData& ssedata, unsigned int oprtr) const override;
            double FindOperatorWeight(SSEData& ssedata, unsigned int leg) const override;
            double FindOperatorWeight(int oprtr_type) const override;
            bool IsNegativeOperator(SSEData& ssedata, unsigned int leg) const override;
            void SpecialUpdates(SSEData& ssedata) override;
            bool IsDiagonal(SSEData& ssedata, unsigned int leg) const override;
            bool IsDiagonal(int oprtr_type) const override;
            Spin RandomSpin() override;
        protected:
            bool IsNegativeOperator(std::vector<Spin> spn_state) const override;
            int FindOperatorType(std::vector<Spin> spn_state) const override;
            std::vector<Spin> Propagate(std::vector<Spin> spn_state, unsigned int oprtr_type) const override;
        private:
            unsigned int length; //!< Number of sites (and bonds) of the lattice.
            std::uniform_int_distribution<int> rdm_spin_sz; //!< Generates a random spin z-component (-1 or 1).

    };

    //! Implementation of a linear one dimensionnal lattice structure.
    class LatticeDiamond : public Lattice{
        public:
            LatticeDiamond(unsigned int length, double j2j1);
            virtual ~LatticeDiamond(){};
            std::vector<unsigned int> GetSpinsIndex(unsigned int num_bond) const override;
            unsigned int GetSpinNumber() const override;
            unsigned int GetBondNumber() const override;
            unsigned int HowManyOperators() const override;
            bool IsLoopStart(SSEData& ssedata, unsigned int leg) const override;
            Spin GetLoopType(SSEData& ssedata, unsigned int leg) override;
            std::vector<double> FindPathProbs(SSEData& ssedata, unsigned int leg, const Spin& loop_type) const override;
            int FindOperatorType(SSEData& ssedata, unsigned int leg) const override;
            int FindBondType(SSEData& ssedata, unsigned int bond) const override;
            std::vector<Spin> PropagateBond(SSEData& ssedata, unsigned int oprtr) const override;
            double FindOperatorWeight(SSEData& ssedata, unsigned int leg) const override;
            double FindOperatorWeight(int oprtr_type) const override;
            bool IsNegativeOperator(SSEData& ssedata, unsigned int leg) const override;
            void SpecialUpdates(SSEData& ssedata) override;
            void ASiteUpdate(SSEData& ssedata);
            void BSiteUpdate(SSEData& ssedata);
            bool IsDiagonal(SSEData& ssedata, unsigned int leg) const override;
            bool IsDiagonal(int oprtr_type) const override;
            Spin RandomSpin() override;
        protected:
            bool IsNegativeOperator(std::vector<Spin> spn_state) const override;
            int FindOperatorType(std::vector<Spin> spn_state) const override;
            std::vector<Spin> Propagate(std::vector<Spin> spn_state, unsigned int oprtr_type) const override;
        private:
            unsigned int length; //!< Number of sites (and bonds) of the lattice.
            double epsilon; //!< Stores the ratio of the coupling constants J1 and J2.
            std::uniform_int_distribution<unsigned int> rdm_spin_norm; //!< Provides a random spin norm (0,1) for the qmc::LatticeDiamond::RandomSpin function.
            std::uniform_int_distribution<int> rdm_spin_sz; //!< Provides a random spin z-component (-1,0,1) to the qmc::LatticeDiamond::RandomSpin function.
    };
}

#endif // QMC_LATTICE_HPP_

