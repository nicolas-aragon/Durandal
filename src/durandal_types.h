#ifndef DURANDAL_TYPES_H
#define DURANDAL_TYPES_H

#include "ffi_vec.h"
#include <vector>
#include <NTL/mat_GF2.h>
#include <NTL/mat_GF2E.h>

using namespace std;
using namespace NTL;

typedef struct publicKey
{
	ffi_vec H1;
	ffi_vec H2;
	vector<ffi_vec> T;
	vector<ffi_vec> Tprime;
} publicKey;

typedef struct secretKey
{
	ffi_vec E;
	vector<ffi_vec> S;
	vector<ffi_vec> Sprime;
} secretKey;

typedef struct offSignature
{
	ffi_vec W;
	ffi_vec F;
	ffi_vec EF;

	ffi_vec y1;
	ffi_vec y2;
	ffi_vec x;
	
	ffi_vec U;
	mat_GF2 D;

	mat_GF2E Smat;
} offSignature;

typedef struct signature
{
	ffi_vec z;
	ffi_vec F;
	unsigned char c[C_BYTES];
	unsigned char p[P_BYTES];
} signature;

#endif
