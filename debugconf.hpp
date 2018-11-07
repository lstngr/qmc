#ifndef QMC_DEBUG_ENGINE_HPP_
#define QMC_DEBUG_ENGINE_HPP_

#define LOG_LEVEL 1
/*
 * Choose a log-level
 *-1 -> None
 * 0 -> Error
 * 1 -> Warning
 * 2 -> Info
 * 3 -> More info
 * 4 -> MOAR INFO
 */

#define LOG(lvl,x) if(lvl<=LOG_LEVEL && lvl>-1){ std::cerr << x; }

#endif //QMC_DEBUG_ENGINE_HPP_
