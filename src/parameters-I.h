#ifndef PARAMETERS_H
#define PARAMETERS_H

#define PARAM_M 241
#define PARAM_N 202
#define PARAM_K 101

#define PARAM_L 4
#define PARAM_LPRIME 1

#define PARAM_D 6
#define PARAM_R 6

#define PARAM_W 57
#define PARAM_LAMBDA 12

#define PARAM_Q 2 //Cannot be modified in this implementation

#define VEC_D_BYTES 181
#define VEC_K_BYTES 3043
#define VEC_N_BYTES 6086

#define SHA512_BYTES 64
//Number of bytes to represent c and p in a basis of F
#define C_BYTES 76
#define P_BYTES 303

//Polynomial P defining the ideal
#define NMODCOEFFS 5
#define MODCOEFFS {0, 1, 6, 7, 101}

#endif
