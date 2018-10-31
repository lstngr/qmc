#include "omanager.hpp"
#include <iomanip>
#include "mcstep.hpp"

using namespace qmc;

OManager::OManager(short unsigned int stream_type, const MCStep* sim_engine, std::string filename): sim_engine_(sim_engine), stream_type_(stream_type) {
    switch(stream_type){
        case 0:
            out_buffer_ = std::cout.rdbuf();
            break;
        case 1:
            out_streamfile_ = std::ofstream();
            out_streamfile_.open(filename);
            out_buffer_ = out_streamfile_.rdbuf();
            break;
        default:
            out_buffer_ = std::cout.rdbuf();
    }

    out_stream_ = new std::ostream(out_buffer_);
}

OManager::~OManager() {
    *out_stream_ << std::endl << std::endl;
    if(out_streamfile_.is_open()){
        out_streamfile_.close();
    }
    delete out_stream_;
}

void OManager::PrintObservables(unsigned int sweep_index) {
    const std::vector<Observable*> obsrv(sim_engine_->GetObservables());
    *out_stream_ << std::setw(10) << std::left << sweep_index << " "/*<< std::setw(7) << sim_engine_->GetMeasureBins()*/;
    for(auto it=obsrv.begin(); it!=obsrv.end(); ++it) {
        *out_stream_ << std::setw(15) << std::setprecision(14) << std::left << (*it)->Measure() << " ";
    }
    *out_stream_ << std::endl;
}

void OManager::PrintInitLine() {
    const std::vector<Observable*> obsrv(sim_engine_->GetObservables());
    *out_stream_ << std::setw(10) << std::left << "# Sweeps";
    for(auto it=obsrv.begin(); it!=obsrv.end(); ++it) {
        *out_stream_ << std::setw(15) << std::left << (*it)->GetName() << " ";
    }
    *out_stream_ << std::endl;
}
