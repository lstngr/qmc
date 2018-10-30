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

#define LOG(lvl,x) if(lvl<=LOG_LEVEL){ std::cerr << x; }
