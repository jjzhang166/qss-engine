#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "pkg_file.h"
#include <common/utils.h>


#include <common/model.h>
#include <qss/qss_model.h>
#include <classic/classic_model.h>

double __dT;
double __ax;
double __ay;
double __r;
double __dx;
double __dy;

void
MOD_settings(SD_simulationSettings settings)
{
	 settings->debug = 0;
	 settings->parallel = FALSE;
	 settings->hybrid = FALSE;
	 settings->method = 4;
}

void
MOD_definition(int i, double *x, double *d, double *alg, double t, double *dx)
{
	int j = 0;
	switch(i)
	{
		case 0:
			dx[1] = (-x[0]*__ax/__dx)+(-x[0]*__ay/__dy)+__r*(pow(x[0],2.0))*(1.0-x[0]);
			return;
		case 1000:
			dx[1] = (-x[3000]*__ax/__dx)+(-x[3000]+x[0])*__ay/__dy+__r*(pow(x[3000],2.0))*(2.0-x[3000]);
			return;
		case 2000:
			dx[1] = (-x[6000]*__ax/__dx)+(-x[6000]+x[3000])*__ay/__dy+__r*(pow(x[6000],2.0))*(2.0-x[6000]);
			return;
		case 3000:
			dx[1] = (-x[9000]*__ax/__dx)+(-x[9000]+x[6000])*__ay/__dy+__r*(pow(x[9000],2.0))*(2.0-x[9000]);
			return;
		case 4000:
			dx[1] = (-x[12000]*__ax/__dx)+(-x[12000]+x[9000])*__ay/__dy+__r*(pow(x[12000],2.0))*(2.0-x[12000]);
			return;
		case 5000:
			dx[1] = (-x[15000]*__ax/__dx)+(-x[15000]+x[12000])*__ay/__dy+__r*(pow(x[15000],2.0))*(2.0-x[15000]);
			return;
		case 6000:
			dx[1] = (-x[18000]*__ax/__dx)+(-x[18000]+x[15000])*__ay/__dy+__r*(pow(x[18000],2.0))*(2.0-x[18000]);
			return;
		case 7000:
			dx[1] = (-x[21000]*__ax/__dx)+(-x[21000]+x[18000])*__ay/__dy+__r*(pow(x[21000],2.0))*(2.0-x[21000]);
			return;
		case 8000:
			dx[1] = (-x[24000]*__ax/__dx)+(-x[24000]+x[21000])*__ay/__dy+__r*(pow(x[24000],2.0))*(2.0-x[24000]);
			return;
		case 9000:
			dx[1] = (-x[27000]*__ax/__dx)+(-x[27000]+x[24000])*__ay/__dy+__r*(pow(x[27000],2.0))*(2.0-x[27000]);
			return;
		default:
			j = i;
			if(j >=1 && j <= 999)
			{
				dx[1] = (-x[(j) * 3]+x[(j-1) * 3])*__ax/__dx+(-x[(j) * 3]*__ay/__dy)+__r*(pow(x[(j) * 3],2.0))*(1.0-x[(j) * 3]);
			}
			j = i-1000;
			if(j >=1 && j <= 999)
			{
				dx[1] = (-x[(j+1000) * 3]+x[(j+999) * 3])*__ax/__dx+(-x[(j+1000) * 3]+x[(j) * 3])*__ay/__dy+__r*(pow(x[(j+1000) * 3],2.0))*(1.0-x[(j+1000) * 3]);
			}
			j = i-2000;
			if(j >=1 && j <= 999)
			{
				dx[1] = (-x[(j+2000) * 3]+x[(j+1999) * 3])*__ax/__dx+(-x[(j+2000) * 3]+x[(j+1000) * 3])*__ay/__dy+__r*(pow(x[(j+2000) * 3],2.0))*(1.0-x[(j+2000) * 3]);
			}
			j = i-3000;
			if(j >=1 && j <= 999)
			{
				dx[1] = (-x[(j+3000) * 3]+x[(j+2999) * 3])*__ax/__dx+(-x[(j+3000) * 3]+x[(j+2000) * 3])*__ay/__dy+__r*(pow(x[(j+3000) * 3],2.0))*(1.0-x[(j+3000) * 3]);
			}
			j = i-4000;
			if(j >=1 && j <= 999)
			{
				dx[1] = (-x[(j+4000) * 3]+x[(j+3999) * 3])*__ax/__dx+(-x[(j+4000) * 3]+x[(j+3000) * 3])*__ay/__dy+__r*(pow(x[(j+4000) * 3],2.0))*(1.0-x[(j+4000) * 3]);
			}
			j = i-5000;
			if(j >=1 && j <= 999)
			{
				dx[1] = (-x[(j+5000) * 3]+x[(j+4999) * 3])*__ax/__dx+(-x[(j+5000) * 3]+x[(j+4000) * 3])*__ay/__dy+__r*(pow(x[(j+5000) * 3],2.0))*(1.0-x[(j+5000) * 3]);
			}
			j = i-6000;
			if(j >=1 && j <= 999)
			{
				dx[1] = (-x[(j+6000) * 3]+x[(j+5999) * 3])*__ax/__dx+(-x[(j+6000) * 3]+x[(j+5000) * 3])*__ay/__dy+__r*(pow(x[(j+6000) * 3],2.0))*(1.0-x[(j+6000) * 3]);
			}
			j = i-7000;
			if(j >=1 && j <= 999)
			{
				dx[1] = (-x[(j+7000) * 3]+x[(j+6999) * 3])*__ax/__dx+(-x[(j+7000) * 3]+x[(j+6000) * 3])*__ay/__dy+__r*(pow(x[(j+7000) * 3],2.0))*(1.0-x[(j+7000) * 3]);
			}
			j = i-8000;
			if(j >=1 && j <= 999)
			{
				dx[1] = (-x[(j+8000) * 3]+x[(j+7999) * 3])*__ax/__dx+(-x[(j+8000) * 3]+x[(j+7000) * 3])*__ay/__dy+__r*(pow(x[(j+8000) * 3],2.0))*(1.0-x[(j+8000) * 3]);
			}
			j = i-9000;
			if(j >=1 && j <= 999)
			{
				dx[1] = (-x[(j+9000) * 3]+x[(j+8999) * 3])*__ax/__dx+(-x[(j+9000) * 3]+x[(j+8000) * 3])*__ay/__dy+__r*(pow(x[(j+9000) * 3],2.0))*(1.0-x[(j+9000) * 3]);
			}
	}
}

void
MOD_dependencies(int i, double *x, double *d, double *alg, double t, double *der, int *map)
{
	int j = 0;
	switch(i)
	{
		case 0:
			der[0 + 1] = (-x[0]*__ax/__dx)+(-x[0]*__ay/__dy)+__r*(pow(x[0],2.0))*(1.0-x[0]);
			der[3000 + 1] = (-x[3000]*__ax/__dx)+(-x[3000]+x[0])*__ay/__dy+__r*(pow(x[3000],2.0))*(2.0-x[3000]);
			break;
		case 1000:
			der[3000 + 1] = (-x[3000]*__ax/__dx)+(-x[3000]+x[0])*__ay/__dy+__r*(pow(x[3000],2.0))*(2.0-x[3000]);
			der[6000 + 1] = (-x[6000]*__ax/__dx)+(-x[6000]+x[3000])*__ay/__dy+__r*(pow(x[6000],2.0))*(2.0-x[6000]);
			break;
		case 2000:
			der[6000 + 1] = (-x[6000]*__ax/__dx)+(-x[6000]+x[3000])*__ay/__dy+__r*(pow(x[6000],2.0))*(2.0-x[6000]);
			der[9000 + 1] = (-x[9000]*__ax/__dx)+(-x[9000]+x[6000])*__ay/__dy+__r*(pow(x[9000],2.0))*(2.0-x[9000]);
			break;
		case 3000:
			der[9000 + 1] = (-x[9000]*__ax/__dx)+(-x[9000]+x[6000])*__ay/__dy+__r*(pow(x[9000],2.0))*(2.0-x[9000]);
			der[12000 + 1] = (-x[12000]*__ax/__dx)+(-x[12000]+x[9000])*__ay/__dy+__r*(pow(x[12000],2.0))*(2.0-x[12000]);
			break;
		case 4000:
			der[12000 + 1] = (-x[12000]*__ax/__dx)+(-x[12000]+x[9000])*__ay/__dy+__r*(pow(x[12000],2.0))*(2.0-x[12000]);
			der[15000 + 1] = (-x[15000]*__ax/__dx)+(-x[15000]+x[12000])*__ay/__dy+__r*(pow(x[15000],2.0))*(2.0-x[15000]);
			break;
		case 5000:
			der[15000 + 1] = (-x[15000]*__ax/__dx)+(-x[15000]+x[12000])*__ay/__dy+__r*(pow(x[15000],2.0))*(2.0-x[15000]);
			der[18000 + 1] = (-x[18000]*__ax/__dx)+(-x[18000]+x[15000])*__ay/__dy+__r*(pow(x[18000],2.0))*(2.0-x[18000]);
			break;
		case 6000:
			der[18000 + 1] = (-x[18000]*__ax/__dx)+(-x[18000]+x[15000])*__ay/__dy+__r*(pow(x[18000],2.0))*(2.0-x[18000]);
			der[21000 + 1] = (-x[21000]*__ax/__dx)+(-x[21000]+x[18000])*__ay/__dy+__r*(pow(x[21000],2.0))*(2.0-x[21000]);
			break;
		case 7000:
			der[21000 + 1] = (-x[21000]*__ax/__dx)+(-x[21000]+x[18000])*__ay/__dy+__r*(pow(x[21000],2.0))*(2.0-x[21000]);
			der[24000 + 1] = (-x[24000]*__ax/__dx)+(-x[24000]+x[21000])*__ay/__dy+__r*(pow(x[24000],2.0))*(2.0-x[24000]);
			break;
		case 8000:
			der[24000 + 1] = (-x[24000]*__ax/__dx)+(-x[24000]+x[21000])*__ay/__dy+__r*(pow(x[24000],2.0))*(2.0-x[24000]);
			der[27000 + 1] = (-x[27000]*__ax/__dx)+(-x[27000]+x[24000])*__ay/__dy+__r*(pow(x[27000],2.0))*(2.0-x[27000]);
			break;
		case 9000:
			der[27000 + 1] = (-x[27000]*__ax/__dx)+(-x[27000]+x[24000])*__ay/__dy+__r*(pow(x[27000],2.0))*(2.0-x[27000]);
			break;
	}
	j = i+1;
	if(j >=1 && j <= 999)
	{
		der[(j) * 3 + 1] = (-x[(j) * 3]+x[(j-1) * 3])*__ax/__dx+(-x[(j) * 3]*__ay/__dy)+__r*(pow(x[(j) * 3],2.0))*(1.0-x[(j) * 3]);
	}
	j = i;
	if(j >=1 && j <= 999)
	{
		der[(j) * 3 + 1] = (-x[(j) * 3]+x[(j-1) * 3])*__ax/__dx+(-x[(j) * 3]*__ay/__dy)+__r*(pow(x[(j) * 3],2.0))*(1.0-x[(j) * 3]);
		der[(j+1000) * 3 + 1] = (-x[(j+1000) * 3]+x[(j+999) * 3])*__ax/__dx+(-x[(j+1000) * 3]+x[(j) * 3])*__ay/__dy+__r*(pow(x[(j+1000) * 3],2.0))*(1.0-x[(j+1000) * 3]);
	}
	j = i-999;
	if(j >=1 && j <= 999)
	{
		der[(j+1000) * 3 + 1] = (-x[(j+1000) * 3]+x[(j+999) * 3])*__ax/__dx+(-x[(j+1000) * 3]+x[(j) * 3])*__ay/__dy+__r*(pow(x[(j+1000) * 3],2.0))*(1.0-x[(j+1000) * 3]);
	}
	j = i-1000;
	if(j >=1 && j <= 999)
	{
		der[(j+1000) * 3 + 1] = (-x[(j+1000) * 3]+x[(j+999) * 3])*__ax/__dx+(-x[(j+1000) * 3]+x[(j) * 3])*__ay/__dy+__r*(pow(x[(j+1000) * 3],2.0))*(1.0-x[(j+1000) * 3]);
		der[(j+2000) * 3 + 1] = (-x[(j+2000) * 3]+x[(j+1999) * 3])*__ax/__dx+(-x[(j+2000) * 3]+x[(j+1000) * 3])*__ay/__dy+__r*(pow(x[(j+2000) * 3],2.0))*(1.0-x[(j+2000) * 3]);
	}
	j = i-1999;
	if(j >=1 && j <= 999)
	{
		der[(j+2000) * 3 + 1] = (-x[(j+2000) * 3]+x[(j+1999) * 3])*__ax/__dx+(-x[(j+2000) * 3]+x[(j+1000) * 3])*__ay/__dy+__r*(pow(x[(j+2000) * 3],2.0))*(1.0-x[(j+2000) * 3]);
	}
	j = i-2000;
	if(j >=1 && j <= 999)
	{
		der[(j+2000) * 3 + 1] = (-x[(j+2000) * 3]+x[(j+1999) * 3])*__ax/__dx+(-x[(j+2000) * 3]+x[(j+1000) * 3])*__ay/__dy+__r*(pow(x[(j+2000) * 3],2.0))*(1.0-x[(j+2000) * 3]);
		der[(j+3000) * 3 + 1] = (-x[(j+3000) * 3]+x[(j+2999) * 3])*__ax/__dx+(-x[(j+3000) * 3]+x[(j+2000) * 3])*__ay/__dy+__r*(pow(x[(j+3000) * 3],2.0))*(1.0-x[(j+3000) * 3]);
	}
	j = i-2999;
	if(j >=1 && j <= 999)
	{
		der[(j+3000) * 3 + 1] = (-x[(j+3000) * 3]+x[(j+2999) * 3])*__ax/__dx+(-x[(j+3000) * 3]+x[(j+2000) * 3])*__ay/__dy+__r*(pow(x[(j+3000) * 3],2.0))*(1.0-x[(j+3000) * 3]);
	}
	j = i-3000;
	if(j >=1 && j <= 999)
	{
		der[(j+3000) * 3 + 1] = (-x[(j+3000) * 3]+x[(j+2999) * 3])*__ax/__dx+(-x[(j+3000) * 3]+x[(j+2000) * 3])*__ay/__dy+__r*(pow(x[(j+3000) * 3],2.0))*(1.0-x[(j+3000) * 3]);
		der[(j+4000) * 3 + 1] = (-x[(j+4000) * 3]+x[(j+3999) * 3])*__ax/__dx+(-x[(j+4000) * 3]+x[(j+3000) * 3])*__ay/__dy+__r*(pow(x[(j+4000) * 3],2.0))*(1.0-x[(j+4000) * 3]);
	}
	j = i-3999;
	if(j >=1 && j <= 999)
	{
		der[(j+4000) * 3 + 1] = (-x[(j+4000) * 3]+x[(j+3999) * 3])*__ax/__dx+(-x[(j+4000) * 3]+x[(j+3000) * 3])*__ay/__dy+__r*(pow(x[(j+4000) * 3],2.0))*(1.0-x[(j+4000) * 3]);
	}
	j = i-4000;
	if(j >=1 && j <= 999)
	{
		der[(j+4000) * 3 + 1] = (-x[(j+4000) * 3]+x[(j+3999) * 3])*__ax/__dx+(-x[(j+4000) * 3]+x[(j+3000) * 3])*__ay/__dy+__r*(pow(x[(j+4000) * 3],2.0))*(1.0-x[(j+4000) * 3]);
		der[(j+5000) * 3 + 1] = (-x[(j+5000) * 3]+x[(j+4999) * 3])*__ax/__dx+(-x[(j+5000) * 3]+x[(j+4000) * 3])*__ay/__dy+__r*(pow(x[(j+5000) * 3],2.0))*(1.0-x[(j+5000) * 3]);
	}
	j = i-4999;
	if(j >=1 && j <= 999)
	{
		der[(j+5000) * 3 + 1] = (-x[(j+5000) * 3]+x[(j+4999) * 3])*__ax/__dx+(-x[(j+5000) * 3]+x[(j+4000) * 3])*__ay/__dy+__r*(pow(x[(j+5000) * 3],2.0))*(1.0-x[(j+5000) * 3]);
	}
	j = i-5000;
	if(j >=1 && j <= 999)
	{
		der[(j+5000) * 3 + 1] = (-x[(j+5000) * 3]+x[(j+4999) * 3])*__ax/__dx+(-x[(j+5000) * 3]+x[(j+4000) * 3])*__ay/__dy+__r*(pow(x[(j+5000) * 3],2.0))*(1.0-x[(j+5000) * 3]);
		der[(j+6000) * 3 + 1] = (-x[(j+6000) * 3]+x[(j+5999) * 3])*__ax/__dx+(-x[(j+6000) * 3]+x[(j+5000) * 3])*__ay/__dy+__r*(pow(x[(j+6000) * 3],2.0))*(1.0-x[(j+6000) * 3]);
	}
	j = i-5999;
	if(j >=1 && j <= 999)
	{
		der[(j+6000) * 3 + 1] = (-x[(j+6000) * 3]+x[(j+5999) * 3])*__ax/__dx+(-x[(j+6000) * 3]+x[(j+5000) * 3])*__ay/__dy+__r*(pow(x[(j+6000) * 3],2.0))*(1.0-x[(j+6000) * 3]);
	}
	j = i-6000;
	if(j >=1 && j <= 999)
	{
		der[(j+6000) * 3 + 1] = (-x[(j+6000) * 3]+x[(j+5999) * 3])*__ax/__dx+(-x[(j+6000) * 3]+x[(j+5000) * 3])*__ay/__dy+__r*(pow(x[(j+6000) * 3],2.0))*(1.0-x[(j+6000) * 3]);
		der[(j+7000) * 3 + 1] = (-x[(j+7000) * 3]+x[(j+6999) * 3])*__ax/__dx+(-x[(j+7000) * 3]+x[(j+6000) * 3])*__ay/__dy+__r*(pow(x[(j+7000) * 3],2.0))*(1.0-x[(j+7000) * 3]);
	}
	j = i-6999;
	if(j >=1 && j <= 999)
	{
		der[(j+7000) * 3 + 1] = (-x[(j+7000) * 3]+x[(j+6999) * 3])*__ax/__dx+(-x[(j+7000) * 3]+x[(j+6000) * 3])*__ay/__dy+__r*(pow(x[(j+7000) * 3],2.0))*(1.0-x[(j+7000) * 3]);
	}
	j = i-7000;
	if(j >=1 && j <= 999)
	{
		der[(j+7000) * 3 + 1] = (-x[(j+7000) * 3]+x[(j+6999) * 3])*__ax/__dx+(-x[(j+7000) * 3]+x[(j+6000) * 3])*__ay/__dy+__r*(pow(x[(j+7000) * 3],2.0))*(1.0-x[(j+7000) * 3]);
		der[(j+8000) * 3 + 1] = (-x[(j+8000) * 3]+x[(j+7999) * 3])*__ax/__dx+(-x[(j+8000) * 3]+x[(j+7000) * 3])*__ay/__dy+__r*(pow(x[(j+8000) * 3],2.0))*(1.0-x[(j+8000) * 3]);
	}
	j = i-7999;
	if(j >=1 && j <= 999)
	{
		der[(j+8000) * 3 + 1] = (-x[(j+8000) * 3]+x[(j+7999) * 3])*__ax/__dx+(-x[(j+8000) * 3]+x[(j+7000) * 3])*__ay/__dy+__r*(pow(x[(j+8000) * 3],2.0))*(1.0-x[(j+8000) * 3]);
	}
	j = i-8000;
	if(j >=1 && j <= 999)
	{
		der[(j+8000) * 3 + 1] = (-x[(j+8000) * 3]+x[(j+7999) * 3])*__ax/__dx+(-x[(j+8000) * 3]+x[(j+7000) * 3])*__ay/__dy+__r*(pow(x[(j+8000) * 3],2.0))*(1.0-x[(j+8000) * 3]);
		der[(j+9000) * 3 + 1] = (-x[(j+9000) * 3]+x[(j+8999) * 3])*__ax/__dx+(-x[(j+9000) * 3]+x[(j+8000) * 3])*__ay/__dy+__r*(pow(x[(j+9000) * 3],2.0))*(1.0-x[(j+9000) * 3]);
	}
	j = i-8999;
	if(j >=1 && j <= 999)
	{
		der[(j+9000) * 3 + 1] = (-x[(j+9000) * 3]+x[(j+8999) * 3])*__ax/__dx+(-x[(j+9000) * 3]+x[(j+8000) * 3])*__ay/__dy+__r*(pow(x[(j+9000) * 3],2.0))*(1.0-x[(j+9000) * 3]);
	}
	j = i-9000;
	if(j >=1 && j <= 999)
	{
		der[(j+9000) * 3 + 1] = (-x[(j+9000) * 3]+x[(j+8999) * 3])*__ax/__dx+(-x[(j+9000) * 3]+x[(j+8000) * 3])*__ay/__dy+__r*(pow(x[(j+9000) * 3],2.0))*(1.0-x[(j+9000) * 3]);
	}
}

void
MOD_output(int i, double *x, double *d, double *alg, double t, double *out)
{
	int j = 0;
	j = i;
	if(j >=0 && j <= 4)
	{
		out[0] = x[(200*j+4000) * 3];
	}
}

void
QSS_initializeDataStructs(QSS_simulator simulator)
{
	int *discretes = (int*)malloc(2*sizeof(int));
	int *outputs = (int*)malloc(5*sizeof(int));
	int *states = (int*)malloc(10000*sizeof(int));
	int i0;
	int i;
	int j = 0;
	simulator->data = QSS_Data(10000,2,0,0,0,"ar2d2");
QSS_data modelData = simulator->data;

	// Allocate main data structures.
	__dT = 1.000000000000000020816682e-03;
	__ax = 1.0;
	__ay = 1.0;
	__r = 1000.0;
	__dx = 10.0/1000;
	__dy = 10.0/10;
	for(i = 0; i <= 999;i++)
	{
		modelData->x[i * 3] = 0.0;
	}
	for(i = 0; i <= 999;i++)
	{
		modelData->x[(i+9000) * 3] = 0.0;
	}
	for(i = 0; i <= 999;i++)
	{
		modelData->x[(i+1000) * 3] = 0.0;
	}
	for(i = 0; i <= 999;i++)
	{
		modelData->x[(i+2000) * 3] = 0.0;
	}
	for(i = 0; i <= 999;i++)
	{
		modelData->x[(i+3000) * 3] = 0.0;
	}
	for(i = 0; i <= 999;i++)
	{
		modelData->x[(i+4000) * 3] = 0.0;
	}
	for(i = 0; i <= 999;i++)
	{
		modelData->x[(i+5000) * 3] = 0.0;
	}
	for(i = 0; i <= 999;i++)
	{
		modelData->x[(i+6000) * 3] = 0.0;
	}
	for(i = 0; i <= 999;i++)
	{
		modelData->x[(i+7000) * 3] = 0.0;
	}
	for(i = 0; i <= 999;i++)
	{
		modelData->x[(i+8000) * 3] = 0.0;
	}
	// Initialize model code.
	for(i0 = 0; i0 <= 99; i0++)
	{
		modelData->x[(i0) * 3] = 1.0;
	}
	modelData->nDS[0] = 1;
	modelData->nDS[1000] = 2;
	modelData->nDS[2000] = 2;
	modelData->nDS[3000] = 2;
	modelData->nDS[4000] = 2;
	modelData->nDS[5000] = 2;
	modelData->nDS[6000] = 2;
	modelData->nDS[7000] = 2;
	modelData->nDS[8000] = 2;
	modelData->nDS[9000] = 2;
	for(i = 1; i <= 999; i++)
	{
		modelData->nDS[i] = 2;
		modelData->nDS[i+1000] = 3;
		modelData->nDS[i+2000] = 3;
		modelData->nDS[i+3000] = 3;
		modelData->nDS[i+4000] = 3;
		modelData->nDS[i+5000] = 3;
		modelData->nDS[i+6000] = 3;
		modelData->nDS[i+7000] = 3;
		modelData->nDS[i+8000] = 3;
		modelData->nDS[i+9000] = 3;
	}
	modelData->nSD[0]++;
	modelData->nSD[0]++;
	modelData->nSD[1000]++;
	modelData->nSD[1000]++;
	modelData->nSD[2000]++;
	modelData->nSD[2000]++;
	modelData->nSD[3000]++;
	modelData->nSD[3000]++;
	modelData->nSD[4000]++;
	modelData->nSD[4000]++;
	modelData->nSD[5000]++;
	modelData->nSD[5000]++;
	modelData->nSD[6000]++;
	modelData->nSD[6000]++;
	modelData->nSD[7000]++;
	modelData->nSD[7000]++;
	modelData->nSD[8000]++;
	modelData->nSD[8000]++;
	modelData->nSD[9000]++;
	for(i = 1; i <= 999; i++)
	{
		modelData->nSD[i-1]++;
		modelData->nSD[i]++;
		modelData->nSD[i]++;
		modelData->nSD[i+999]++;
		modelData->nSD[i+1000]++;
		modelData->nSD[i+1000]++;
		modelData->nSD[i+1999]++;
		modelData->nSD[i+2000]++;
		modelData->nSD[i+2000]++;
		modelData->nSD[i+2999]++;
		modelData->nSD[i+3000]++;
		modelData->nSD[i+3000]++;
		modelData->nSD[i+3999]++;
		modelData->nSD[i+4000]++;
		modelData->nSD[i+4000]++;
		modelData->nSD[i+4999]++;
		modelData->nSD[i+5000]++;
		modelData->nSD[i+5000]++;
		modelData->nSD[i+5999]++;
		modelData->nSD[i+6000]++;
		modelData->nSD[i+6000]++;
		modelData->nSD[i+6999]++;
		modelData->nSD[i+7000]++;
		modelData->nSD[i+7000]++;
		modelData->nSD[i+7999]++;
		modelData->nSD[i+8000]++;
		modelData->nSD[i+8000]++;
		modelData->nSD[i+8999]++;
		modelData->nSD[i+9000]++;
	}
	QSS_allocDataMatrix(modelData);
	// Initialize model data.
	// Initialize model time.
		cleanVector(states,0,10000);

	modelData->DS[0][states[0]++] = 0;
	modelData->DS[1000][states[1000]++] = 0;
	modelData->DS[1000][states[1000]++] = 1000;
	modelData->DS[2000][states[2000]++] = 1000;
	modelData->DS[2000][states[2000]++] = 2000;
	modelData->DS[3000][states[3000]++] = 2000;
	modelData->DS[3000][states[3000]++] = 3000;
	modelData->DS[4000][states[4000]++] = 3000;
	modelData->DS[4000][states[4000]++] = 4000;
	modelData->DS[5000][states[5000]++] = 4000;
	modelData->DS[5000][states[5000]++] = 5000;
	modelData->DS[6000][states[6000]++] = 5000;
	modelData->DS[6000][states[6000]++] = 6000;
	modelData->DS[7000][states[7000]++] = 6000;
	modelData->DS[7000][states[7000]++] = 7000;
	modelData->DS[8000][states[8000]++] = 7000;
	modelData->DS[8000][states[8000]++] = 8000;
	modelData->DS[9000][states[9000]++] = 8000;
	modelData->DS[9000][states[9000]++] = 9000;
	for(i = 1; i <= 999; i++)
	{
		modelData->DS[i][states[i]++] = i-1;
		modelData->DS[i][states[i]++] = i;
		modelData->DS[i+1000][states[i+1000]++] = i;
		modelData->DS[i+1000][states[i+1000]++] = i+999;
		modelData->DS[i+1000][states[i+1000]++] = i+1000;
		modelData->DS[i+2000][states[i+2000]++] = i+1000;
		modelData->DS[i+2000][states[i+2000]++] = i+1999;
		modelData->DS[i+2000][states[i+2000]++] = i+2000;
		modelData->DS[i+3000][states[i+3000]++] = i+2000;
		modelData->DS[i+3000][states[i+3000]++] = i+2999;
		modelData->DS[i+3000][states[i+3000]++] = i+3000;
		modelData->DS[i+4000][states[i+4000]++] = i+3000;
		modelData->DS[i+4000][states[i+4000]++] = i+3999;
		modelData->DS[i+4000][states[i+4000]++] = i+4000;
		modelData->DS[i+5000][states[i+5000]++] = i+4000;
		modelData->DS[i+5000][states[i+5000]++] = i+4999;
		modelData->DS[i+5000][states[i+5000]++] = i+5000;
		modelData->DS[i+6000][states[i+6000]++] = i+5000;
		modelData->DS[i+6000][states[i+6000]++] = i+5999;
		modelData->DS[i+6000][states[i+6000]++] = i+6000;
		modelData->DS[i+7000][states[i+7000]++] = i+6000;
		modelData->DS[i+7000][states[i+7000]++] = i+6999;
		modelData->DS[i+7000][states[i+7000]++] = i+7000;
		modelData->DS[i+8000][states[i+8000]++] = i+7000;
		modelData->DS[i+8000][states[i+8000]++] = i+7999;
		modelData->DS[i+8000][states[i+8000]++] = i+8000;
		modelData->DS[i+9000][states[i+9000]++] = i+8000;
		modelData->DS[i+9000][states[i+9000]++] = i+8999;
		modelData->DS[i+9000][states[i+9000]++] = i+9000;
	}
		cleanVector(states,0,10000);

	modelData->SD[0][states[0]++] = 0;
	modelData->SD[0][states[0]++] = 1000;
	modelData->SD[1000][states[1000]++] = 1000;
	modelData->SD[1000][states[1000]++] = 2000;
	modelData->SD[2000][states[2000]++] = 2000;
	modelData->SD[2000][states[2000]++] = 3000;
	modelData->SD[3000][states[3000]++] = 3000;
	modelData->SD[3000][states[3000]++] = 4000;
	modelData->SD[4000][states[4000]++] = 4000;
	modelData->SD[4000][states[4000]++] = 5000;
	modelData->SD[5000][states[5000]++] = 5000;
	modelData->SD[5000][states[5000]++] = 6000;
	modelData->SD[6000][states[6000]++] = 6000;
	modelData->SD[6000][states[6000]++] = 7000;
	modelData->SD[7000][states[7000]++] = 7000;
	modelData->SD[7000][states[7000]++] = 8000;
	modelData->SD[8000][states[8000]++] = 8000;
	modelData->SD[8000][states[8000]++] = 9000;
	modelData->SD[9000][states[9000]++] = 9000;
	for(i = 1; i <= 999; i++)
	{
		modelData->SD[i-1][states[i-1]++] = i;
		modelData->SD[i][states[i]++] = i;
		modelData->SD[i][states[i]++] = i+1000;
		modelData->SD[i+999][states[i+999]++] = i+1000;
		modelData->SD[i+1000][states[i+1000]++] = i+1000;
		modelData->SD[i+1000][states[i+1000]++] = i+2000;
		modelData->SD[i+1999][states[i+1999]++] = i+2000;
		modelData->SD[i+2000][states[i+2000]++] = i+2000;
		modelData->SD[i+2000][states[i+2000]++] = i+3000;
		modelData->SD[i+2999][states[i+2999]++] = i+3000;
		modelData->SD[i+3000][states[i+3000]++] = i+3000;
		modelData->SD[i+3000][states[i+3000]++] = i+4000;
		modelData->SD[i+3999][states[i+3999]++] = i+4000;
		modelData->SD[i+4000][states[i+4000]++] = i+4000;
		modelData->SD[i+4000][states[i+4000]++] = i+5000;
		modelData->SD[i+4999][states[i+4999]++] = i+5000;
		modelData->SD[i+5000][states[i+5000]++] = i+5000;
		modelData->SD[i+5000][states[i+5000]++] = i+6000;
		modelData->SD[i+5999][states[i+5999]++] = i+6000;
		modelData->SD[i+6000][states[i+6000]++] = i+6000;
		modelData->SD[i+6000][states[i+6000]++] = i+7000;
		modelData->SD[i+6999][states[i+6999]++] = i+7000;
		modelData->SD[i+7000][states[i+7000]++] = i+7000;
		modelData->SD[i+7000][states[i+7000]++] = i+8000;
		modelData->SD[i+7999][states[i+7999]++] = i+8000;
		modelData->SD[i+8000][states[i+8000]++] = i+8000;
		modelData->SD[i+8000][states[i+8000]++] = i+9000;
		modelData->SD[i+8999][states[i+8999]++] = i+9000;
		modelData->SD[i+9000][states[i+9000]++] = i+9000;
	}
	simulator->time = QSS_Time(10000,0,0,0,ST_Binary,NULL);

	simulator->output = SD_Output("ar2d2",5,2,10000,NULL,0,0,CI_Step,SD_Memory,MOD_output);
SD_output modelOutput = simulator->output;

	for(i = 0; i <= 4; i++)
	{
		modelOutput->nOS[i] = 1;
		modelOutput->nSO[200*i+4000]++;
	}
	SD_allocOutputMatrix(modelOutput,10000,2);
		cleanVector(states,0,10000);

		cleanVector(outputs,0,5);

	for(i = 0; i <= 4; i++)
	{
		sprintf(modelOutput->variable[i].name,"u5[%d]",200*i+1);
	}
		cleanVector(outputs,0,5);

	for(i = 0; i <= 4; i++)
	{
		modelOutput->SO[200*i+4000][states[200*i+4000]++] = i;
		modelOutput->OS[i][outputs[i]++] = 200*i+4000;
	}
	simulator->model = QSS_Model(MOD_definition,MOD_dependencies,NULL,NULL,NULL);
	free(discretes);
	free(outputs);
	free(states);
}

void
CLC_initializeDataStructs (CLC_simulator simulator)
{
}
