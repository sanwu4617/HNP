#include "gmp_head.h"
mpz_t mpz_pow2[1024];
void gmp_init()
{
    mpz_init_set_str(mpz_pow2[0],"1",10);
    mpz_init_set_str(mpz_pow2[1],"2",10);
    for(int i=2;i<1024;i++)
    {
        mpz_init(mpz_pow2[i]);
        mpz_mul(mpz_pow2[i],mpz_pow2[i-1],mpz_pow2[1]);
    }

}