#include <stdlib.h>
#include <math.h>

#include "../../src/engine/common/random.h"
#include "mmo_math.h"

double mmo_rand(double v)
{
	double y;
	y = v*rand()/RAND_MAX;
	return y;
}

unsigned long mmo_getRandomValue(int n)
{
	return getRandomValue(n);
}

double mmo_exponential(double mu)
{
	return exponential(mu);
}

double mmo_uniform(double a, double b)
{
	return uniform(a,b);
}

double mmo_normal(double sigma)
{
	return normal(sigma);
}

void mmo_shuffle(int *a, int size)
{
	return shuffle(a,size);
}
