#define DEBUG_MODE

#ifdef DEBUG_MODE
#define DEBUG_GEN(x) do { std::cerr << x; } while (0)
#define DEBUG_GEN_BLOCK

#else
#define DEBUG_GEN(x)

#endif //DEBUG_MODE
