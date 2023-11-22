#ifndef NTL_HEAD_H
#define NTL_HEAD_H

#include <NTL/matrix.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/LLL.h>
#include <fstream>
#include <time.h>
NTL_CLIENT
using namespace std;
using namespace NTL;

extern ZZ ZZ_pow2[1024];
extern double d_pow2[1024];
extern RR RR_pow2[1024];

void NTL_init();

#endif