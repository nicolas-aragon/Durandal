#include "signature.h"
#include "durandal_types.h"
#include "ffi_field.h"

#define ITERATIONS 1

int main() {
	//Initializes GF(2^m)
	ffi_field_init();
	//Initializes the polynomial defined the ideal
	ffi_vec_init_mulmod();

	clock_t begin, end;
	double keygenTime = 0;
	double offlineSignatureTime = 0;
	double onlineSignatureTime = 0;
	double verifyTime = 0;

	for(int iter=0 ; iter<ITERATIONS ; iter++) {
		char* message = "Test";
		unsigned int mlen = 4;

		publicKey pk;
		secretKey sk;

		offSignature offlineSignature;
		signature signature;

		begin = clock();
		keygen(pk, sk);
		end = clock();
		keygenTime += (double)(end - begin) / CLOCKS_PER_SEC;

		begin = clock();
		offlineSign(pk, sk, offlineSignature);
		end = clock();
		offlineSignatureTime += (double)(end-begin) / CLOCKS_PER_SEC;

		begin = clock();
		onlineSign(pk, sk, offlineSignature, signature, message, mlen);
		end = clock();
		onlineSignatureTime += (double)(end-begin) / CLOCKS_PER_SEC;

		begin = clock();
		int status = verify(pk, signature, message, mlen);
		end = clock();
		verifyTime += (double)(end-begin) / CLOCKS_PER_SEC;

		if(status == 0) cout << "Valid signature" << endl;
		else cout << "Invalid signature" << endl;
	}

	cout << "Keygen : " << keygenTime / ITERATIONS << endl;
	cout << "Offline signature : " << offlineSignatureTime / ITERATIONS << endl;
	cout << "Online Signature : " << onlineSignatureTime / ITERATIONS << endl;
	cout << "Verification : " << verifyTime / ITERATIONS << endl;
}
