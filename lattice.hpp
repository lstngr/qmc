#ifndef LATTICE_HPP_
#define LATTICE_HPP_

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
     */
    class Lattice{
        public:
            Lattice();
            virtual ~Lattice(){};
            virtual std::vector<unsigned int> GetSpinsIndex(unsigned int num_bond) const=0;
        protected:
        private:

    };

    //! Implementation of a linear one dimensionnal lattice structure.
    /*!
     * Name is pretty self-explanatory
     */
    class Lattice1D : public Lattice{
        public:
            Lattice1D();
            virtual ~Lattice1D(){};
            std::vector<unsigned int> GetSpinsIndex(unsigned int num_bond) const;
        protected:
        private:

    };

}

#endif // LATTICE_HPP_

