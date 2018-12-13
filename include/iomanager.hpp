#ifndef QMC_OMANAGER_HPP_
#define QMC_OMANAGER_HPP_

#include <vector>
#include <ostream>
#include <iostream>
#include <fstream>
#include <map>
// #include "debugconf.hpp"

namespace qmc{

    // Forward declaration for pointers
    class MCStep;

    /*! \brief Static Simulation Information
     *
     * Static (and passive) simulation parameters are stored in this structure.
     * They are initialized by qmc::IManager and are passed to qmc::MCStep.
     */
    struct Configuration{
        Configuration();
        double beta; //!< Inverse temperature of the system
        double expansion_factor; //!< Safety factor for operator strin expansion
        unsigned int measurement_frequency; //!< Number of measures per observable bin
        std::string LoopMethod; //!< Either Wolff or Swendsen-Wang algorithms
        unsigned int num_sweep; //!< Number of sweeps to perform
        int max_sweep; //<! Maximum allowed number of sweeps (prevents runaway due to negative configurations)
        std::string output_dest; //!< Either a std stream or a filename
        std::map<std::string,int> param2num; //!<Maps parameter to an index
    };

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
            void PrintInitLine();
        private:
            std::ostream* out_stream_; //!< Output stream pointing to either std::cout or an output file
            std::streambuf* out_buffer_; //!< Reads a buffer from std::ostream children classes
            std::ofstream out_streamfile_; //!< Output stream to a file.
            const MCStep* sim_engine_; //!< Pointer to the simulation manager to access content to print.
            const short unsigned int stream_type_; //!< Index setting the type of output stream to be used.
    };

    /*! \brief Input manager
     *
     * This class is responsible for reading input files when a simulation is
     * launched.
     */
    class IManager{
        public:
            IManager(std::string filename="default.in");
            virtual ~IManager();
            Configuration GetConfig();
        private:
            void ImportParameter(const std::string& key, const std::string& value);
            Configuration sim_conf_; //!< Configuration object storing parameters to be read.
            std::istream* is_file; //!< Input stream pointing to the configuration file.
    };
}

#endif //QMC_OMANAGER_HPP_
