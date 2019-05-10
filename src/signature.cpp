#include "signature.h"
#include "ffi_vec.h"
#include "ffi_elt.h"

#include "hash.h"

#include <stdio.h>
#include <vector>
#include <math.h>

#include <NTL/mat_GF2.h>
#include <NTL/mat_GF2E.h>
#include <NTL/vec_GF2.h>

using namespace std;

//Declaration of additionnal functions
static void get_bit(unsigned char* b, unsigned char c, unsigned int position);
static void set_bit(unsigned char* c, unsigned int position);

bool filtering(const ffi_vec &U, const ffi_vec &E, const ffi_vec &F);
void systemCompute(mat_GF2 &D, const vector<ffi_vec> &S, const ffi_vec &E, const vector<vec_GF2> &mappings, mat_GF2E &Smat);
void unfoldInSupport_precompute(const ffi_vec &support, const int &rank, mat_GF2 &invSystem, vector<int> &columns);
vec_GF2 unfoldInSupport(const ffi_elt &v, const mat_GF2 &invSystem, const int &rank, const vector<int> &columns);

/** Protocol functions **/

void keygen(publicKey &pk, secretKey &sk) {
	//First generate the matrix H, by choosing H1 (left part) and H2 (right part) randomly
	ffi_vec_set_random_full_rank_using_rng(pk.H1, PARAM_K);
	ffi_vec_set_random_full_rank_using_rng(pk.H2, PARAM_K);

	//Generate the support E and the s_i
	ffi_vec_set_random_full_rank_using_rng(sk.E, PARAM_R);

	ffi_vec tmp, tmp1, tmp2;

	//Each s_i is split into two polynomials of length K because of the ideal structure
	for(int i=0 ; i<2*PARAM_L ; i++) {
		ffi_vec_set_random_from_support_using_rng(tmp, PARAM_K, sk.E, PARAM_R);
		sk.S.push_back(tmp);
	}

	for(int i=0 ; i<2*PARAM_LPRIME ; i++) {
		ffi_vec_set_random_from_support_using_rng(tmp, PARAM_K, sk.E, PARAM_R);
		sk.Sprime.push_back(tmp);
	}

	//Compute the corresponding t_i
	for(int i=0 ; i<PARAM_L ; i++) {
		ffi_vec_mul(tmp1, pk.H1, sk.S[2*i], PARAM_K);
		ffi_vec_mul(tmp2, pk.H2, sk.S[2*i+1], PARAM_K);
		ffi_vec_add(tmp, tmp1, tmp2, PARAM_K);
		pk.T.push_back(tmp);
	}

	for(int i=0 ; i<PARAM_LPRIME ; i++) {
		ffi_vec_mul(tmp1, pk.H1, sk.Sprime[2*i], PARAM_K);
		ffi_vec_mul(tmp2, pk.H2, sk.Sprime[2*i+1], PARAM_K);
		ffi_vec_add(tmp, tmp1, tmp2, PARAM_K);
		pk.Tprime.push_back(tmp);
	}
}

void offlineSign(const publicKey &pk, const secretKey &sk, offSignature &offlineSignature) {
	//First randomly choose W and F
	ffi_vec product, w_ef;
	do {
		ffi_vec_set_random_full_rank_using_rng(offlineSignature.W, PARAM_W);
		ffi_vec_set_random_full_rank_using_rng(offlineSignature.F, PARAM_D);

		//Compute the vector space W + EF
		//The order of E and F is important ot get the basis under this form : {f1.e1, f2.e1, ... }
		ffi_vec_tensor_mul(product, sk.E, PARAM_R, offlineSignature.F, PARAM_D);
		ffi_vec_directsum(w_ef, product, PARAM_R * PARAM_D, offlineSignature.W, PARAM_W);
	} while (ffi_vec_get_rank(w_ef, PARAM_W + PARAM_D * PARAM_R) != PARAM_W + PARAM_D * PARAM_R); //Ensure that the vector space is of dimension W + EF

	//Randomly sample y in this vector space
	//Split into two polynomials of length k because of the ideal structure
	ffi_vec_set_random_from_support_using_rng(offlineSignature.y1, PARAM_K, w_ef, PARAM_W + PARAM_D * PARAM_R);
	ffi_vec_set_random_from_support_using_rng(offlineSignature.y2, PARAM_K, w_ef, PARAM_W + PARAM_D * PARAM_R);

	//Compute x
	ffi_vec tmp, tmp2;
	ffi_vec_mul(tmp, pk.H1, offlineSignature.y1, PARAM_K);
	ffi_vec_mul(tmp2, pk.H2, offlineSignature.y2, PARAM_K);
	ffi_vec_add(offlineSignature.x, tmp, tmp2, PARAM_K);

	do {
		//Choose U, a random subspace of EF of dimension rd - lambda, such that no element of the form e.f, e \in E, f \in F, is included in U
		do {
			do {
				ffi_vec_set_random_from_support_using_rng(offlineSignature.U, PARAM_R * PARAM_D - PARAM_LAMBDA, product, PARAM_R * PARAM_D);
			} while (ffi_vec_get_rank(offlineSignature.U, PARAM_R * PARAM_D - PARAM_LAMBDA) != PARAM_R * PARAM_D - PARAM_LAMBDA);
		} while (!filtering(offlineSignature.U, sk.E, offlineSignature.F));

		//Find elements such that they form a complementary basis of U in EF
		//Store this basis - will be used during computation of p
		do {
			ffi_vec_set_random_from_support_using_rng(offlineSignature.EF, PARAM_LAMBDA, product, PARAM_R * PARAM_D);
			offlineSignature.EF.SetLength(PARAM_R * PARAM_D);
			for(int i=0 ; i<PARAM_R * PARAM_D - PARAM_LAMBDA ; i++) {
				offlineSignature.EF[i+PARAM_LAMBDA] = offlineSignature.U[i];
			}
		} while(ffi_vec_get_rank(offlineSignature.EF, PARAM_R * PARAM_D) != PARAM_R * PARAM_D);

		//Unfold each coordinate of product into this new basis to build the change of basis matrix
		//Precomputation for unfolding in EF
		mat_GF2 system;
		vector<int> columns;
		unfoldInSupport_precompute(offlineSignature.EF, PARAM_R * PARAM_D, system, columns);

		mat_GF2 matrix;
		matrix.SetDims(PARAM_R * PARAM_D, PARAM_R * PARAM_D);
		for(int i=0 ; i<PARAM_R * PARAM_D ; i++) {
			vec_GF2 unfoldedProduct = unfoldInSupport(product[i], system, PARAM_R * PARAM_D, columns);
			for(int j=0 ; j<PARAM_R * PARAM_D ; j++) {
				matrix[j][i] = unfoldedProduct[j];
			}
		}

		//Get the lambda first lines to get lines corresponding to the complementary basis of U in EF
		vector<vec_GF2> mappings;
		for(int i=0 ; i<PARAM_LAMBDA; i++) {
			mappings.push_back(matrix[i]);
		}

		//We now need to compute the matrix D, that will be used during the online signature phase to speed up the system solving
		systemCompute(offlineSignature.D, sk.S, sk.E, mappings, offlineSignature.Smat);
	} while(offlineSignature.D.NumRows() == 0); //Repeat until this system has a solution
}

void onlineSign(const publicKey &pk, const secretKey &sk, const offSignature &offlineSignature, signature &signature, char* message, unsigned int messageLen) {
	//Compute c by hashing x, F, and the message mu together
	//We need to obtain PARAM_L * PARAM_K bits to construct c
	
	unsigned char hash[SHA512_BYTES];
	unsigned char hash_input[VEC_D_BYTES + VEC_K_BYTES + messageLen];
	ffi_vec_to_string(hash_input, offlineSignature.F, PARAM_D);
	ffi_vec_to_string(hash_input + VEC_D_BYTES, offlineSignature.x, PARAM_K);
	memcpy(hash_input + VEC_D_BYTES + VEC_K_BYTES, message, messageLen);
	sha512(hash, hash_input, VEC_D_BYTES + VEC_K_BYTES +messageLen);

	memcpy(signature.c, hash, SHA512_BYTES);
	int copiedBytes = SHA512_BYTES;
	while(C_BYTES > copiedBytes) {
		memcpy(hash_input, hash, SHA512_BYTES);
		sha512(hash, hash_input, SHA512_BYTES);

		if(C_BYTES - copiedBytes > SHA512_BYTES) memcpy(signature.c + copiedBytes, hash, SHA512_BYTES);
		else memcpy(signature.c + copiedBytes, hash, C_BYTES - SHA512_BYTES);

		copiedBytes += SHA512_BYTES;
	}

	vec_GF2E cTmp;
	cTmp.SetLength(PARAM_LPRIME * PARAM_K);

	for(int i=0 ; i<PARAM_K * PARAM_LPRIME ; i++) {
		ffi_elt acc;
		ffi_elt_set_zero(acc);
		for(int j=0 ; j<PARAM_D ; j++) {
			unsigned char bit;
			int pos = i * PARAM_D + j;
			get_bit(&bit, signature.c[pos/8], pos % 8);
			if(bit) {
				ffi_elt_add(acc, acc, offlineSignature.F[j]);
			}
		}
		cTmp[i] = acc;
	}

	//Compute z (and p)
	//First compute y + cS'
	//Build the matrix S' explicitly
	mat_GF2E Sprimemat;
	Sprimemat.SetDims(PARAM_LPRIME * PARAM_K, PARAM_N);

	ffi_elt one;
	ffi_elt_set_one(one);

	ffi_vec X;
	ffi_vec_set_zero(X, PARAM_K);
	ffi_vec_set_coeff(X, one, 1);

	for(int i=0 ; i<PARAM_LPRIME ; i++) {
		ffi_vec leftPart, rightPart;
		ffi_vec_set(leftPart, sk.Sprime[2*i], PARAM_K);
		ffi_vec_set(rightPart, sk.Sprime[2*i + 1], PARAM_K);

		for(int k=0 ; k<PARAM_K ; k++) {
				//Copy the left and right part
				for(int pos=0 ; pos<PARAM_K ; pos++) Sprimemat[i*PARAM_K + k][pos] = leftPart[pos];
				for(int pos=0 ; pos<PARAM_K ; pos++) Sprimemat[i*PARAM_K + k][pos + PARAM_K] = rightPart[pos];

			//Compute the ideal shifts
			ffi_vec_mul(leftPart, leftPart, X, PARAM_K);
			ffi_vec_mul(rightPart, rightPart, X, PARAM_K);
		}
	}

	//Use of NTL functions for simplicity
	vec_GF2E zTmp = cTmp * Sprimemat;

	//Add zTmp to y to obtain y + cS'
	ffi_vec_set_zero(signature.z, PARAM_N);

	for(int i=0 ; i<PARAM_K ; i++) {
		ffi_elt acc;

		ffi_elt_set_zero(acc);
		ffi_elt_set(acc, zTmp[i]);
		ffi_elt_add(acc, acc, offlineSignature.y1[i]);
		ffi_vec_set_coeff(signature.z, acc, i);

		ffi_elt_set_zero(acc);
		ffi_elt_set(acc, zTmp[i+PARAM_K]);
		ffi_elt_add(acc, acc, offlineSignature.y2[i]);
		ffi_vec_set_coeff(signature.z, acc, i+PARAM_K);
	}

	//Unfold this vector in a basis of EF to find the coordinates corresponding to the complementary basis
	//The first lambda coordinates correspond to equations when solving for p

	//Precomputation for unfolding in EF
	mat_GF2 system;
	vector<int> columns;
	unfoldInSupport_precompute(offlineSignature.EF, PARAM_R * PARAM_D, system, columns);

	vec_GF2 equations;
	equations.SetLength(PARAM_N * PARAM_LAMBDA);
	for(int i=0 ; i<PARAM_N ; i++) {
		vec_GF2 unfolded = unfoldInSupport(signature.z[i], system, PARAM_R * PARAM_D, columns);
		for(int j=0 ; j<PARAM_LAMBDA ; j++) {
			equations[i*PARAM_LAMBDA + j] = unfolded[j];
		}
	}

	//Compute equations \times D to obtain the coordinates of p
	//Store the coordinates of p in the signature
	vec_GF2 pCoords;
	pCoords = equations * offlineSignature.D;

	//Clear signature.p then copy
	for(int i=0 ; i<P_BYTES ; i++) signature.p[i] = 0;
	for(int i=0 ; i<PARAM_L * PARAM_K * PARAM_D ; i++) {
		if(pCoords[i] == 1) set_bit(signature.p + (i/8), i%8);
	}

	//Build the vector p in Fqm
	vec_GF2E pTmp;
	pTmp.SetLength(PARAM_L * PARAM_K);
	for(int i=0 ; i<PARAM_L*PARAM_K ; i++) {
		ffi_elt acc;
		ffi_elt_set_zero(acc);
		for(int j=0 ; j<PARAM_D ; j++) {
			if(pCoords[i*PARAM_D + j] == 1) {
				ffi_elt_add(acc, acc, offlineSignature.F[j]);
			}
		}
		pTmp[i] = acc;
	}

	vec_GF2E pS;
	pS = pTmp * offlineSignature.Smat;

	for(int i=0 ; i<PARAM_N ; i++) {
		ffi_elt acc;
		ffi_elt_set_zero(acc);
		ffi_elt_add(acc, signature.z[i], pS[i]);
		ffi_vec_set_coeff(signature.z, acc, i);
	}

	//Output (z, F, c, p)
	ffi_vec_set(signature.F, offlineSignature.F, PARAM_D);

	for(int i=0 ; i<C_BYTES ; i++) printf("%.2X", signature.c[i]);
	printf("\n");
}

//Returns 0 if the signature is valid, 1 otherwise
int verify(const publicKey &pk, const signature &signature, char* message, unsigned int messageLen) {
	//First verify the weight of z
	if(ffi_vec_get_rank(signature.z, PARAM_N) > PARAM_W + PARAM_R * PARAM_D - PARAM_LAMBDA) return 1;

	ffi_elt one;
	ffi_elt_set_one(one);

	ffi_vec X;
	ffi_vec_set_zero(X, PARAM_K);
	ffi_vec_set_coeff(X, one, 1);

	//Explicitly construct T and Tprime
	mat_GF2E Tmat;
	Tmat.SetDims(PARAM_K, PARAM_L * PARAM_K);

	for(int i=0 ; i<PARAM_L ; i++) {
		ffi_vec currentShift;
		ffi_vec_set(currentShift, pk.T[i], PARAM_K);

		for(int k=0 ; k<PARAM_K ; k++) {
			for(int pos=0 ; pos<PARAM_K ; pos++) Tmat[pos][i*PARAM_K + k] = currentShift[pos];

			//Compute the ideal shift
			ffi_vec_mul(currentShift, currentShift, X, PARAM_K);
		}
	}

	mat_GF2E Tmatprime;
	Tmatprime.SetDims(PARAM_K, PARAM_LPRIME * PARAM_K);

	for(int i=0 ; i<PARAM_LPRIME ; i++) {
		ffi_vec currentShift;
		ffi_vec_set(currentShift, pk.Tprime[i], PARAM_K);

		for(int k=0 ; k<PARAM_K ; k++) {
			for(int pos=0 ; pos<PARAM_K ; pos++) Tmatprime[pos][i*PARAM_K + k] = currentShift[pos];

			//Compute the ideal shift
			ffi_vec_mul(currentShift, currentShift, X, PARAM_K);
		}
	}

	//Recover c and p in Fqm
	vec_GF2E cTmp, pTmp;
	cTmp.SetLength(PARAM_LPRIME * PARAM_K);
	pTmp.SetLength(PARAM_L * PARAM_K);

	for(int i=0 ; i<PARAM_K * PARAM_LPRIME ; i++) {
		ffi_elt acc;
		ffi_elt_set_zero(acc);
		for(int j=0 ; j<PARAM_D ; j++) {
			unsigned char bit;
			int pos = i * PARAM_D + j;
			get_bit(&bit, signature.c[pos/8], pos % 8);
			if(bit) {
				ffi_elt_add(acc, acc, signature.F[j]);
			}
		}
		cTmp[i] = acc;
	}

	for(int i=0 ; i<PARAM_K * PARAM_L ; i++) {
		ffi_elt acc;
		ffi_elt_set_zero(acc);
		for(int j=0 ; j<PARAM_D ; j++) {
			unsigned char bit;
			int pos = i * PARAM_D + j;
			get_bit(&bit, signature.p[pos/8], pos % 8);
			if(bit) {
				ffi_elt_add(acc, acc, signature.F[j]);
			}
		}
		pTmp[i] = acc;
	}

	//Compute Hz - T'c + Tp
	ffi_vec z1, z2;
	ffi_vec_set_zero(z1, PARAM_K);
	ffi_vec_set_zero(z2, PARAM_K);
	for(int i=0 ; i<PARAM_K ; i++) {
		ffi_vec_set_coeff(z1, signature.z[i], i);
		ffi_vec_set_coeff(z2, signature.z[i+PARAM_K], i);
	}

	//Hz
	ffi_vec Hz, Hz1, Hz2;
	ffi_vec_mul(Hz1, z1, pk.H1, PARAM_K);
	ffi_vec_mul(Hz2, z2, pk.H2, PARAM_K);
	ffi_vec_add(Hz, Hz1, Hz2, PARAM_K);

	vec_GF2E Tprimec;
	Tprimec = Tmatprime * cTmp;

	vec_GF2E Tp;
	Tp = Tmat * pTmp;

	ffi_vec verif;
	ffi_vec_set(verif, Hz, PARAM_K);
	for(int i=0 ; i<PARAM_K ; i++) {
		ffi_elt acc;
		//Only add is implemented since, for q=2, add and sub are equivalent
		ffi_elt_add(acc, verif[i], Tprimec[i]);
		ffi_elt_add(acc, acc, Tp[i]);
		ffi_vec_set_coeff(verif, acc, i);
	}

	unsigned char verifHash[C_BYTES];

	unsigned char hash[SHA512_BYTES];
	unsigned char hash_input[VEC_D_BYTES + VEC_K_BYTES + messageLen];
	ffi_vec_to_string(hash_input, signature.F, PARAM_D);
	ffi_vec_to_string(hash_input + VEC_D_BYTES, verif, PARAM_K);
	memcpy(hash_input + VEC_D_BYTES + VEC_K_BYTES, message, messageLen);
	sha512(hash, hash_input, VEC_D_BYTES + VEC_K_BYTES + messageLen);

	memcpy(verifHash, hash, SHA512_BYTES);
	int copiedBytes = SHA512_BYTES;
	while(C_BYTES > copiedBytes) {
		memcpy(hash_input, hash, SHA512_BYTES);
		sha512(hash, hash_input, SHA512_BYTES);

		if(C_BYTES - copiedBytes > SHA512_BYTES) memcpy(verifHash + copiedBytes, hash, SHA512_BYTES);
		else memcpy(verifHash + copiedBytes, hash, C_BYTES - SHA512_BYTES);

		copiedBytes += SHA512_BYTES;
	}

	for(int i=0 ; i<C_BYTES ; i++) printf("%.2X", verifHash[i]);
	printf("\n");

	//Verify that this hash is equal to c
	if(memcmp(signature.c, verifHash, C_BYTES) == 0) return 0;
	else return 1;
}

















/** Additionnal functions **/

static void get_bit(unsigned char* b, unsigned char c, unsigned int position) {
  *b = (c >> position) & 0x01;
}

static void set_bit(unsigned char* c, unsigned int position) {
  *c = *c | (1 << position);
}

bool filtering(const ffi_vec &U, const ffi_vec &E, const ffi_vec &F) {
	//Enumerate every element of E
	for (int i=0 ; i<pow(2, PARAM_R) ; i++) {
		ffi_elt e;

		for(int j=0 ; j<PARAM_R ; j++) {
			if(i & (1 << j)) ffi_elt_add(e, e, E[j]);
		}

		//Compute e_i.F
		ffi_vec Fi;
		ffi_vec_scalar_mul(Fi, F, e, PARAM_D);

		//Check the intersection with U is null
		ffi_vec inter;
		unsigned int dim;
		ffi_vec_intersection(inter, dim, Fi, PARAM_D, U, PARAM_R * PARAM_D - PARAM_LAMBDA);

		if(dim > 0) return false;
	}

	return true;
}

void systemCompute(mat_GF2 &D, const vector<ffi_vec> &S, const ffi_vec &E, const vector<vec_GF2> &mappings, mat_GF2E &Smat) {
	//We start by writing the matrix lk \times n matrix S including the ideal shifts
	Smat.SetDims(PARAM_L * PARAM_K, PARAM_N);

	ffi_elt one;
	ffi_elt_set_one(one);

	ffi_vec X;
	ffi_vec_set_zero(X, PARAM_K);
	ffi_vec_set_coeff(X, one, 1);

	for(int i=0 ; i<PARAM_L ; i++) {
		ffi_vec leftPart, rightPart;
		ffi_vec_set(leftPart, S[2*i], PARAM_K);
		ffi_vec_set(rightPart, S[2*i + 1], PARAM_K);

		for(int k=0 ; k<PARAM_K ; k++) {
			//Copy the left and right part
			for(int pos=0 ; pos<PARAM_K ; pos++) Smat[i*PARAM_K + k][pos] = leftPart[pos];
			for(int pos=0 ; pos<PARAM_K ; pos++) Smat[i*PARAM_K + k][pos + PARAM_K] = rightPart[pos];

			//Compute the ideal shifts
			ffi_vec_mul(leftPart, leftPart, X, PARAM_K);
			ffi_vec_mul(rightPart, rightPart, X, PARAM_K);
		}
	}

	//Each coordinate in F2m is then unfolded into a d \times rd matrix in F2 as in the formal decoding of LRPC codes
	mat_GF2 fullSystem;
	fullSystem.SetDims(PARAM_L * PARAM_K * PARAM_D, PARAM_N * PARAM_R * PARAM_D);

	//Precomputation for unfolding in E
	mat_GF2 system;
	vector<int> columns;
	unfoldInSupport_precompute(E, PARAM_R, system, columns);

	for(int i=0 ; i<PARAM_K*PARAM_L ; i++) {
		for(int j=0 ; j<PARAM_N ; j++) {
			//Unfold the coordinates
			vec_GF2 coords = unfoldInSupport(Smat[i][j], system, PARAM_R, columns);

			//Write the d \times rd matrix
			for(int k=0 ; k<PARAM_R ; k++) {
				for(int l=0 ; l<PARAM_D ; l++) {
					fullSystem[i*PARAM_D + l][j * PARAM_R * PARAM_D + k * PARAM_D + l] = coords[k];
				}
			}
		}
	}

	//Use the provided mappings to compute the system corresponding to the complementary basis of U in EF
	mat_GF2 finalSystem;
	finalSystem.SetDims(PARAM_L * PARAM_K * PARAM_D, PARAM_N * PARAM_LAMBDA);

	for(int row=0 ; row<PARAM_L * PARAM_K * PARAM_D ; row++) {
		for(int i=0 ; i<PARAM_N ; i++) {
			for(int j=0 ; j<PARAM_LAMBDA ; j++) {
				for(int k=0 ; k<PARAM_R * PARAM_D ; k++) {
					if(mappings[j][k] == 1) {
						//Add the corresponding column
						finalSystem[row][i*PARAM_LAMBDA + j] += fullSystem[row][i*PARAM_R*PARAM_D + k];
					}
				}
			}
		}
	}

	//Finally invert this system to compute the matrix D
	GF2 det;
	inv(det, D, finalSystem);
}


//Precomputes the inverse of the system to unfold a vector in its support
void unfoldInSupport_precompute(const ffi_vec &support, const int &rank, mat_GF2 &invSystem, vector<int> &columns) {
	//Build the system
	mat_GF2 system;
	system.SetDims(rank, rank);

	unsigned int offset;

	//We need to find an invertible square matrix in the rank \times m matrix (support unfolded in F2^m)
	//We store the chosen columns in an vector
	int currentCol = 0;

	for(int column=0 ; column<PARAM_M ; column++) {
		//Copy column into currentCol
		for(int j=0 ; j<rank ; j++) {
		   system[j][currentCol] = coeff(rep(support[j]), column);
		}

		mat_GF2 testRank = system;

		if(gauss(testRank) == currentCol + 1) {
			columns.push_back(column);
			currentCol++;
		}

		if(currentCol == rank) break;
	}

	//Inverse this system
	GF2 det;
	inv(det, invSystem, system);
}

//Unfolds an element of F2m into coordinates in the given basis of its support
//Builds the system each time : can be improved a lot
vec_GF2 unfoldInSupport(const ffi_elt &v, const mat_GF2 &invSystem, const int &rank, const vector<int> &columns) {
	//Solve the system
	vec_GF2 equations;
	equations.SetLength(rank);

	for(int j=0 ; j<rank ; j++) {
		equations[j] = coeff(rep(v), columns[j]);
	}

	//This function makes the assumption that the system has a solution
	vec_GF2 solutions;
	solutions = equations * invSystem;

  return solutions;
}