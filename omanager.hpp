#ifndef QMC_OMANAGER_HPP_
#define QMC_OMANAGER_HPP_

#include <vector>
#include <ostream>
#include <iostream>
#include <fstream>
// #include "debugconf.hpp"

namespace qmc{

    class MCStep;

    /*! \brief Manager class for output streams
     *
     * This class enables writing simulation data to either standard output
     * streams (like cout or cerr), or to specific files. It is meant to manage
     * all sources of information from the simulation, that is, mainly:
     * - Values resulting of physical observables measurements
     * - States (operators/spins) of the system
     * - Debbuging information
     *
     * Common formatting options are built-in (line- or columnwise).
     */
    class OManager{
        public:
            OManager(short unsigned int stream_type, const MCStep* sim_engine, std::string filename="data_out.out");
            virtual ~OManager();
            void PrintObservables(unsigned int sweep_index=0);
            void PrintSpins();
            void PrintOperators();
        private:
            std::ostream* out_stream_;
            std::streambuf* out_buffer_;
            std::ofstream out_streamfile_;
            const MCStep* sim_engine_;
            const short unsigned int stream_type_;
            bool out_stream_open_;
            // TODO STRUCT
            // std::string obsrv_filename_;
            // std::string oprtr_filename_;
            // std::string spins_filename_;
            // bool show_obsrv_;
            // bool show_oprtr_;
            // bool show_spins_;
            // END TODO STRUCT
    };
}

#endif //QMC_OMANAGER_HPP_
