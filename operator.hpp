#ifndef QMC_OPERATOR_HPP_
#define QMC_OPERATOR_HPP_

#include "ssedata.hpp"

namespace qmc{

    class Operator{
        public:
            Operator(){};
            virtual ~Operator(){};
            virtual double Measure(SSEData* const sim_data) const = 0;
        protected:
        private:
    };

};

#endif //QMC_OPERATOR_HPP_
