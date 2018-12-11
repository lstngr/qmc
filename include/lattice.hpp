#ifndef QMC_LATTICE_HPP_
#define QMC_LATTICE_HPP_

#include <vector>

namespace qmc{

    class SSEData;

    struct Spin{
        Spin(unsigned int tsn, int tsz):two_norm(tsn),two_sz(tsz){}
        unsigned int two_norm;
        int two_sz;
    };

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

            virtual bool IsLoopStart(SSEData& ssedata, unsigned int leg) const; // Def always true
            virtual int  GetLoopType(SSEData& ssedata, unsigned int leg) const = 0;
            virtual std::vector<double> GetPathProbs(SSEData& ssedata, unsigned int leg) const = 0;
            virtual unsigned int GetOperatorType(SSEData& ssedata, unsigned int leg) const = 0;
            virtual double GetOperatorWeight(SSEData& ssedata, unsigned int leg) const = 0;
            virtual bool IsNegativeOperator(SSEData& ssedata, unsigned int leg) const = 0;
            virtual void SpecialUpdates(SSEData& ssedata) const = 0;

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
            bool IsLoopStart(SSEData& ssedata, unsigned int leg) const override; // Def always true
            int  GetLoopType(SSEData& ssedata, unsigned int leg) const override;
            std::vector<double> GetPathProbs(SSEData& ssedata, unsigned int leg) const override;
            unsigned int GetOperatorType(SSEData& ssedata, unsigned int leg) const override;
            double GetOperatorWeight(SSEData& ssedata, unsigned int leg) const override;
            bool IsNegativeOperator(SSEData& ssedata, unsigned int leg) const override;
            void SpecialUpdates(SSEData& ssedata) const override;
        protected:
        private:
            unsigned int length; //!< Number of sites (and bonds) of the lattice.
    };

}

#endif // QMC_LATTICE_HPP_

