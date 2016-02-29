#ifndef MMO_MATH_H
#define MMO_MATH_H

double mmo_rand(double v);

unsigned long mmo_getRandomValue(int n);

double mmo_exponential(double mu);

double mmo_uniform(double a, double b);

double mmo_normal(double sigma);

void mmo_shuffle(int *a, int size);

#endif // MMO_MATH_H
