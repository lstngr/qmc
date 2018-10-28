#ifndef QMC_OBSERVABLE_HPP_
#define QMC_OBSERVABLE_HPP_

#include "ssedata.hpp"

namespace qmc{

    class Observable{
        public:
            Observable(){};
            virtual ~Observable(){};
            virtual void UpdateObservable(SSEData* const sim_data)=0;
            virtual double Measure() const = 0;
    };

    class ObservableHamilton : public Observable{
        public:
            ObservableHamilton(double beta);
            virtual ~ObservableHamilton(){};
            virtual void UpdateObservable(SSEData* const sim_data);
            virtual double Measure() const;
        protected:
            std::vector<size_t> num_oprtr_; /*!<Stores the number of inserted operators for every iteration.*/
            double beta_;
    };

};

#endif //QMC_OBSERVABLE_HPP_
