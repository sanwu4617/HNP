#ifndef GMP_HEAD_H
#define GMP_HEAD_H

#include <fplll.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <chrono>
using namespace std;
using namespace fplll;
using namespace std::chrono;

extern mpz_t mpz_pow2[1024];

void gmp_init();

#endif