#include "NTL_head.h"
ZZ ZZ_pow2[1024];
RR RR_pow2[1024];
double d_pow2[1024];
void NTL_init()
{
	ZZ_pow2[0] = to_ZZ("1");
	for (int i = 1; i < 1024; i++)
	{
		ZZ_pow2[i] = 2 * ZZ_pow2[i - 1];
		RR_pow2[i] = MakeRR(ZZ_pow2[i],0);
	}
	d_pow2[0] = 1;
	for (int i = 1; i < 1024; i++)
	{
		d_pow2[i] = 2 * d_pow2[i - 1];
	}
}