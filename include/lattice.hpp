#ifndef QMC_LATTICE_HPP_
#define QMC_LATTICE_HPP_

#include <vector>
#include <ostream>

namespace qmc{

    class SSEData;

    struct Spin{
        Spin(int tsn, int tsz):two_norm(tsn),two_sz(tsz){}
        void Raise(int amount=1){two_sz+=2*amount;}
        void Lower(int amount=1){Raise(-amount);}
        void Promote(int amount=1){two_norm+=2*amount;}
        void Demote (int amount=1){Promote(-amount);}
        int two_norm;
        int two_sz;
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
            virtual unsigned int HowManyOperators() const = 0;
            virtual unsigned int MaxTwoS() const = 0;
            virtual bool IsLoopStart(SSEData& ssedata, unsigned int leg) const = 0;
            virtual int  GetLoopType(SSEData& ssedata, unsigned int leg) const;
            virtual std::vector<double> FindPathProbs(SSEData& ssedata, unsigned int leg) const = 0;
            virtual int FindOperatorType(SSEData& ssedata, unsigned int leg) const = 0;
            virtual int FindOperatorType(std::vector<Spin> spn_state) const = 0;
            virtual double FindOperatorWeight(SSEData& ssedata, unsigned int leg) const = 0;
            virtual double FindOperatorWeight(int oprtr_type) const = 0;
            virtual bool IsNegativeOperator(SSEData& ssedata, unsigned int leg) const = 0;
            virtual void SpecialUpdates(SSEData& ssedata) const = 0;
            virtual bool IsDiagonal(SSEData& ssedata, unsigned int leg) const = 0;
            virtual bool IsDiagonal(int oprtr_type) const = 0;
            virtual std::vector<Spin> Propagate(std::vector<Spin> spn_state, unsigned int oprtr_type) const = 0;

        protected:
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
            unsigned int MaxTwoS() const override{return 1;}
            bool IsLoopStart(SSEData& ssedata, unsigned int leg) const override;
            std::vector<double> FindPathProbs(SSEData& ssedata, unsigned int leg) const override;
            int FindOperatorType(SSEData& ssedata, unsigned int leg) const override;
            int FindOperatorType(std::vector<Spin> spn_state) const override;
            double FindOperatorWeight(SSEData& ssedata, unsigned int leg) const override;
            double FindOperatorWeight(int oprtr_type) const override;
            bool IsNegativeOperator(SSEData& ssedata, unsigned int leg) const override;
            void SpecialUpdates(SSEData& ssedata) const override;
            bool IsDiagonal(SSEData& ssedata, unsigned int leg) const override;
            bool IsDiagonal(int oprtr_type) const override;
            std::vector<Spin> Propagate(std::vector<Spin> spn_state, unsigned int oprtr_type) const override;
        protected:
        private:
            unsigned int length; //!< Number of sites (and bonds) of the lattice.
    };

}

#endif // QMC_LATTICE_HPP_

