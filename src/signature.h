#ifndef SIGNATURE_H
#define SIGNATURE_H

#include "durandal_types.h"

void keygen(publicKey &pk, secretKey &sk);
void offlineSign(const publicKey &pk, const secretKey &sk, offSignature &offlineSignature);
void onlineSign(const publicKey &pk, const secretKey &sk, const offSignature &offlineSignature, signature &signature, char* message, unsigned int messageLen);
int verify(const publicKey &pk, const signature &signature, char* message, unsigned int messageLen);

#endif