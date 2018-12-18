#include "iomanager.hpp"
#include <iomanip>
#include <sstream>
#include "mcstep.hpp"

using namespace qmc;

/* \brief Configuration constructor
 *
 * Build a default configuration with reasonable values. Also builds a table
 * matching parameters and index for ease-of-use in
 * \ref qmc::IManager "import procedures".
 */
Configuration::Configuration()
    :beta(10),expansion_factor(0.2),measurement_frequency(10),
    LoopMethod("wolff"),num_sweep(1000),output_dest("data.out")
{
    param2num = std::map<std::string,int>();
    param2num.insert(std::pair<std::string,int> ("beta",0));
    param2num.insert(std::pair<std::string,int> ("expansion_factor",1));
    param2num.insert(std::pair<std::string,int> ("measurement_frequency",2));
    param2num.insert(std::pair<std::string,int> ("LoopMethod",3));
    param2num.insert(std::pair<std::string,int> ("num_sweep",4));
    param2num.insert(std::pair<std::string,int> ("output_dest",5));
}

/*! \brief Output Manager constructor
 *
 * Builds an OManager object. It mainly manages streams creation and
 * configuration. To remain general in the fashion the data is output, we do not
 * define a precise out-stream type, but use a pointer to std::ostream.
 * \param stream_type Index of the stream type to use.
 * \param sim_engine Pointer to a simulation object to fetch data from.
 * \param filename String of the filename to write to.
 */
OManager::OManager(short unsigned int stream_type, const MCStep* sim_engine, std::string filename): sim_engine_(sim_engine), stream_type_(stream_type) {
    switch(stream_type){
        case 0:
            out_buffer_ = std::cout.rdbuf();
            break;
        case 1:
            out_streamfile_ = std::ofstream();
            out_streamfile_.open(filename, std::ios_base::app);
            out_buffer_ = out_streamfile_.rdbuf();
            break;
        default:
            out_buffer_ = std::cout.rdbuf();
    }

    out_stream_ = new std::ostream(out_buffer_);
}

//! Closes and destroys output streams if necessary
OManager::~OManager() {
    *out_stream_ << std::endl << std::endl;
    if(out_streamfile_.is_open()){
        out_streamfile_.close();
    }
    delete out_stream_;
}

/*! \brief Prints observables from sim_engine_
 *
 * Fetches a reference to all observables emulated by qmc::MCStep and displays
 * them in a column format.
 * \param sweep_index Gotten from sim_engine_, gives current sweep iteration.
 */
void OManager::PrintObservables(unsigned int sweep_index) {
    // TODO -- Harmonize. Either call sim_engine_ to get sweep_idx or remove
    // sim_engine all together.
    const std::vector<Observable*> obsrv(sim_engine_->GetObservables());
    *out_stream_ << std::setw(10) << std::right << sweep_index << ";"/*<< std::setw(7) << sim_engine_->GetMeasureBins()*/;
    for(auto it=obsrv.begin(); it!=obsrv.end(); ++it) {
        *out_stream_ << std::setw(15) << std::setprecision(10) << std::right << (*it)->Measure() << ";";
    }
    *out_stream_ << std::endl;
}

//! Prints an initialization line, usually before the first sweeps
void OManager::PrintInitLine() {
    const std::vector<Observable*> obsrv(sim_engine_->GetObservables());
    *out_stream_ << std::left << "#" << std::setw(10) << std::right << "Sweeps;";
    for(auto it=obsrv.begin(); it!=obsrv.end(); ++it) {
        *out_stream_ << std::setw(15) << std::right << (*it)->GetName() << ";";
    }
    *out_stream_ << std::endl;
}

/*! \brief IManager constructor
 *
 *  Opens an input stream to the specified configuration
 *  \param filename Name of the file to write to.
 */
IManager::IManager(std::string filename) {
    is_file = new std::ifstream(filename);
}

//! Destroys the input stream
IManager::~IManager() {
    delete is_file;
}

/*! \brief Reads an input file and imports its data
 *
 * This function reads the previously set input file line by line and imports
 * the input values in the sim_conf_ object. It overwrites the default values
 * that were set by the class constructor.
 * \return A configuration structure with the updated parameters.
 */
Configuration IManager::GetConfig() {
    std::string read_line;
    while(std::getline(*is_file,read_line)) {
        std::istringstream is_line(read_line);
        std::string key;
        if(std::getline(is_line,key,'=')) {
            std::string value;
            if(std::getline(is_line,value)) {
                ImportParameter(key,value);
            }
        }
    }
    return sim_conf_;
}

/*! \brief Given a key and value, set it in sim_conf_.
 *
 * Using two strings gotten from the \ref qmc::IManager::GetConfig "parser",
 * calls the conversion map stored in sim_conf_ to associate the key string to
 * its integer index. From here, import the value (with correct type conversion)
 * to sim_conf.
 * \param key Key's name
 * \param value String value of the key
 */
void IManager::ImportParameter(const std::string& key, const std::string& value) {
    int key_idx = sim_conf_.param2num.find(key)->second;
    switch(key_idx) {
        case 0:
            sim_conf_.beta = std::atof(value.c_str());
            break;
        case 1:
            sim_conf_.expansion_factor = std::atof(value.c_str());
            break;
        case 2:
            sim_conf_.measurement_frequency = std::stoi(value);
            break;
        case 3:
            sim_conf_.LoopMethod = value;
            break;
        case 4:
            sim_conf_.num_sweep = std::stoi(value);
            break;
        case 5:
            sim_conf_.output_dest = value;
            break;
    }
}
