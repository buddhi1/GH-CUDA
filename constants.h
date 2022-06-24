#define EPSILON_POSITIONS 9
#define EPSILON  0.000000001 
#define DEBUG_MODE 1
#define DEBUG_TIMING 1

// Smaller datasets
// #define xThreadPerBlock 4
// #define yThreadPerBlock 1
// #define yBlockPerGrid 1
// #define MAX_POLY2_SIZE 4
// #define SHARED_MEMORY_PADDING 1

// intermediate datasets
// #define xThreadPerBlock 16
// #define yThreadPerBlock 1
// #define yBlockPerGrid 4
// #define MAX_POLY2_SIZE 512
// #define SHARED_MEMORY_PADDING 1

// larger datasets
#define xThreadPerBlock 32
#define yThreadPerBlock 1
#define yBlockPerGrid 48
#define MAX_POLY2_SIZE 1024
#define SHARED_MEMORY_PADDING 1