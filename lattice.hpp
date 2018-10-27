#ifndef QMC_LATTICE_HPP_
#define QMC_LATTICE_HPP_

#include <vector>

namespace qmc{

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
             */
            virtual std::vector<unsigned int> GetSpinsIndex(unsigned int num_bond) const = 0;

            //! Returns the total number of spins
            /*!
             * Returns the number of spins (or sites) enclosed in the Lattice
             * structure. This function is adapted polymorphically according to
             * the lattice type.
             */
            virtual unsigned int GetSpinNumber() const = 0;

            //! Returns the total number of bonds
            /*!
             * Returns the number of bonds (linking sites) enclosed in the
             * Lattice structure. This function is adapted polymorphically
             * according to the lattice type.
             */
            virtual unsigned int GetBondNumber() const = 0;
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
            std::vector<unsigned int> GetSpinsIndex(unsigned int num_bond) const;
            unsigned int GetSpinNumber() const;
            unsigned int GetBondNumber() const;
        protected:
        private:
            unsigned int length; //!< Number of sites (and bonds) of the lattice.
    };

}

#endif // QMC_LATTICE_HPP_

