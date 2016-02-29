#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <common/utils.h>


#include <common/model.h>
#include <qss/qss_model.h>
#include <classic/classic_model.h>

double __K0;
double __rho0;
double __c;
double __Z0;
double __u0;
double __beta;
double __x0;
double __gam;

void
MOD_settings(SD_simulationSettings settings)
{
	 settings->debug = 0;
	 settings->parallel = FALSE;
	 settings->hybrid = FALSE;
	 settings->method = 5;
}

void
MOD_definition(int i, double *x, double *d, double *alg, double t, double *dx)
{
	int j = 0;
	switch(i)
	{
		case 0:
			dx[1] = __K0*((-x[800])+x[1596])*200;
			return;
		case 200:
			dx[1] = 1.0/__rho0*((-x[4])+x[0])*200;
			return;
		case 199:
			dx[1] = __K0*((-x[1596])+x[1592])*200;
			return;
		case 399:
			dx[1] = 1.0/__rho0*((-x[0])+x[796])*200;
			return;
		case 400:
			dx[1] = (__u0-__c)*((-x[1604])+x[1600])*200;
			return;
		case 600:
			dx[1] = (__u0+__c)*((-x[2400])+x[3196])*200;
			return;
		case 599:
			dx[1] = (__u0-__c)*((-x[1600])+x[2396])*200;
			return;
		case 799:
			dx[1] = (__u0+__c)*((-x[3196])+x[3192])*200;
			return;
		default:
			j = i;
			if(j >=1 && j <= 198)
			{
				dx[1] = __K0*((-x[(j+200) * 4])+x[(j+199) * 4])*200;
			}
			j = i-200;
			if(j >=1 && j <= 198)
			{
				dx[1] = 1.0/__rho0*((-x[(j+1) * 4])+x[(j) * 4])*200;
			}
			j = i-400;
			if(j >=1 && j <= 198)
			{
				dx[1] = (__u0-__c)*((-x[(j+401) * 4])+x[(j+400) * 4])*200;
			}
			j = i-600;
			if(j >=1 && j <= 198)
			{
				dx[1] = (__u0+__c)*((-x[(j+600) * 4])+x[(j+599) * 4])*200;
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
			der[800 + 1] = 1.0/__rho0*((-x[4])+x[0])*200;
			der[1596 + 1] = 1.0/__rho0*((-x[0])+x[796])*200;
			return;
		case 1:
			der[800 + 1] = 1.0/__rho0*((-x[4])+x[0])*200;
			break;
		case 199:
			der[1596 + 1] = 1.0/__rho0*((-x[0])+x[796])*200;
			break;
		case 200:
			der[0 + 1] = __K0*((-x[800])+x[1596])*200;
			break;
		case 398:
			der[796 + 1] = __K0*((-x[1596])+x[1592])*200;
			break;
		case 399:
			der[0 + 1] = __K0*((-x[800])+x[1596])*200;
			der[796 + 1] = __K0*((-x[1596])+x[1592])*200;
			return;
		case 400:
			der[1600 + 1] = (__u0-__c)*((-x[1604])+x[1600])*200;
			der[2396 + 1] = (__u0-__c)*((-x[1600])+x[2396])*200;
			return;
		case 401:
			der[1600 + 1] = (__u0-__c)*((-x[1604])+x[1600])*200;
			break;
		case 599:
			der[2396 + 1] = (__u0-__c)*((-x[1600])+x[2396])*200;
			break;
		case 600:
			der[2400 + 1] = (__u0+__c)*((-x[2400])+x[3196])*200;
			break;
		case 798:
			der[3196 + 1] = (__u0+__c)*((-x[3196])+x[3192])*200;
			break;
		case 799:
			der[2400 + 1] = (__u0+__c)*((-x[2400])+x[3196])*200;
			der[3196 + 1] = (__u0+__c)*((-x[3196])+x[3192])*200;
			return;
	}
	j = i;
	if(j >=1 && j <= 198)
	{
		der[(j+200) * 4 + 1] = 1.0/__rho0*((-x[(j+1) * 4])+x[(j) * 4])*200;
	}
	j = i-1;
	if(j >=1 && j <= 198)
	{
		der[(j+200) * 4 + 1] = 1.0/__rho0*((-x[(j+1) * 4])+x[(j) * 4])*200;
	}
	j = i-199;
	if(j >=1 && j <= 198)
	{
		der[(j) * 4 + 1] = __K0*((-x[(j+200) * 4])+x[(j+199) * 4])*200;
	}
	j = i-200;
	if(j >=1 && j <= 198)
	{
		der[(j) * 4 + 1] = __K0*((-x[(j+200) * 4])+x[(j+199) * 4])*200;
	}
	j = i-400;
	if(j >=1 && j <= 198)
	{
		der[(j+400) * 4 + 1] = (__u0-__c)*((-x[(j+401) * 4])+x[(j+400) * 4])*200;
	}
	j = i-401;
	if(j >=1 && j <= 198)
	{
		der[(j+400) * 4 + 1] = (__u0-__c)*((-x[(j+401) * 4])+x[(j+400) * 4])*200;
	}
	j = i-599;
	if(j >=1 && j <= 198)
	{
		der[(j+600) * 4 + 1] = (__u0+__c)*((-x[(j+600) * 4])+x[(j+599) * 4])*200;
	}
	j = i-600;
	if(j >=1 && j <= 198)
	{
		der[(j+600) * 4 + 1] = (__u0+__c)*((-x[(j+600) * 4])+x[(j+599) * 4])*200;
	}
}

void
MOD_output(int i, double *x, double *d, double *alg, double t, double *out)
{
	int j = 0;
	j = i;
	if(j >=0 && j <= 4)
	{
		out[0] = x[(40*j) * 4];
	}
}

void
QSS_initializeDataStructs(QSS_simulator simulator)
{
	int *outputs = (int*)malloc(5*sizeof(int));
	int *states = (int*)malloc(800*sizeof(int));
	int i0;
	int i;
	int j = 0;
	simulator->data = QSS_Data(800,0,0,0,0,"acousticsRiemDes");
QSS_data modelData = simulator->data;

	// Allocate main data structures.
	__K0 = 9.000000000000000000000000e+00;
	__rho0 = 1.000000000000000000000000e+00;
	__c = 3.000000000000000000000000e+00;
	__Z0 = 3.000000000000000000000000e+00;
	__u0 = 0.000000000000000000000000e+00;
	__beta = 1000.0;
	__x0 = 5.000000000000000000000000e-01;
	__gam = 0.000000000000000000000000e+00;
	// Initialize model code.
	for(i0 = 0; i0 <= 199; i0++)
	{
		modelData->x[(i0) * 4] = pow(2.717999999999999971578291e+00,(-__beta*pow((1.000000000000000000000000e+00*(i0+1)/200-__x0),2.0)))*cos(__gam*(1.000000000000000000000000e+00*(i0+1)/200-__x0));
		modelData->x[(i0+200) * 4] = 0.0;
		modelData->x[(i0+400) * 4] = 1.0/2.0/__Z0*((-modelData->x[(i0) * 4])+modelData->x[(i0+200) * 4]);
		modelData->x[(i0+600) * 4] = 1.0/2.0/__Z0*(modelData->x[(i0) * 4]+modelData->x[(i0+200) * 4]);
	}
	modelData->nDS[0] = 2;
	modelData->nDS[200] = 2;
	for(i = 1; i <= 198; i++)
	{
		modelData->nDS[i] = 2;
		modelData->nDS[i+200] = 2;
	}
	modelData->nDS[199] = 2;
	modelData->nDS[399] = 2;
	modelData->nDS[400] = 2;
	modelData->nDS[600] = 2;
	for(i = 1; i <= 198; i++)
	{
		modelData->nDS[i+400] = 2;
		modelData->nDS[i+600] = 2;
	}
	modelData->nDS[599] = 2;
	modelData->nDS[799] = 2;
	modelData->nSD[200]++;
	modelData->nSD[399]++;
	modelData->nSD[0]++;
	modelData->nSD[1]++;
	for(i = 1; i <= 198; i++)
	{
		modelData->nSD[i+199]++;
		modelData->nSD[i+200]++;
		modelData->nSD[i]++;
		modelData->nSD[i+1]++;
	}
	modelData->nSD[398]++;
	modelData->nSD[399]++;
	modelData->nSD[0]++;
	modelData->nSD[199]++;
	modelData->nSD[400]++;
	modelData->nSD[401]++;
	modelData->nSD[600]++;
	modelData->nSD[799]++;
	for(i = 1; i <= 198; i++)
	{
		modelData->nSD[i+400]++;
		modelData->nSD[i+401]++;
		modelData->nSD[i+599]++;
		modelData->nSD[i+600]++;
	}
	modelData->nSD[400]++;
	modelData->nSD[599]++;
	modelData->nSD[798]++;
	modelData->nSD[799]++;
	QSS_allocDataMatrix(modelData);
	// Initialize model data.
	// Initialize model time.
		cleanVector(states,0,800);

	modelData->DS[0][states[0]++] = 200;
	modelData->DS[0][states[0]++] = 399;
	modelData->DS[200][states[200]++] = 0;
	modelData->DS[200][states[200]++] = 1;
	for(i = 1; i <= 198; i++)
	{
		modelData->DS[i][states[i]++] = i+199;
		modelData->DS[i][states[i]++] = i+200;
		modelData->DS[i+200][states[i+200]++] = i;
		modelData->DS[i+200][states[i+200]++] = i+1;
	}
	modelData->DS[199][states[199]++] = 398;
	modelData->DS[199][states[199]++] = 399;
	modelData->DS[399][states[399]++] = 0;
	modelData->DS[399][states[399]++] = 199;
	modelData->DS[400][states[400]++] = 400;
	modelData->DS[400][states[400]++] = 401;
	modelData->DS[600][states[600]++] = 600;
	modelData->DS[600][states[600]++] = 799;
	for(i = 1; i <= 198; i++)
	{
		modelData->DS[i+400][states[i+400]++] = i+400;
		modelData->DS[i+400][states[i+400]++] = i+401;
		modelData->DS[i+600][states[i+600]++] = i+599;
		modelData->DS[i+600][states[i+600]++] = i+600;
	}
	modelData->DS[599][states[599]++] = 400;
	modelData->DS[599][states[599]++] = 599;
	modelData->DS[799][states[799]++] = 798;
	modelData->DS[799][states[799]++] = 799;
		cleanVector(states,0,800);

	modelData->SD[200][states[200]++] = 0;
	modelData->SD[399][states[399]++] = 0;
	modelData->SD[0][states[0]++] = 200;
	modelData->SD[1][states[1]++] = 200;
	for(i = 1; i <= 198; i++)
	{
		modelData->SD[i+199][states[i+199]++] = i;
		modelData->SD[i+200][states[i+200]++] = i;
		modelData->SD[i][states[i]++] = i+200;
		modelData->SD[i+1][states[i+1]++] = i+200;
	}
	modelData->SD[398][states[398]++] = 199;
	modelData->SD[399][states[399]++] = 199;
	modelData->SD[0][states[0]++] = 399;
	modelData->SD[199][states[199]++] = 399;
	modelData->SD[400][states[400]++] = 400;
	modelData->SD[401][states[401]++] = 400;
	modelData->SD[600][states[600]++] = 600;
	modelData->SD[799][states[799]++] = 600;
	for(i = 1; i <= 198; i++)
	{
		modelData->SD[i+400][states[i+400]++] = i+400;
		modelData->SD[i+401][states[i+401]++] = i+400;
		modelData->SD[i+599][states[i+599]++] = i+600;
		modelData->SD[i+600][states[i+600]++] = i+600;
	}
	modelData->SD[400][states[400]++] = 599;
	modelData->SD[599][states[599]++] = 599;
	modelData->SD[798][states[798]++] = 799;
	modelData->SD[799][states[799]++] = 799;
	simulator->time = QSS_Time(800,0,0,0,ST_Binary,NULL);

	simulator->output = SD_Output("acousticsRiemDes",5,0,800,NULL,0,0,CI_Step,SD_Memory,MOD_output);
SD_output modelOutput = simulator->output;

	for(i = 0; i <= 4; i++)
	{
		modelOutput->nOS[i] = 1;
		modelOutput->nSO[40*i]++;
	}
	SD_allocOutputMatrix(modelOutput,800,0);
		cleanVector(states,0,800);

		cleanVector(outputs,0,5);

	for(i = 0; i <= 4; i++)
	{
		sprintf(modelOutput->variable[i].name,"p[%d]",40*i+1);
	}
		cleanVector(outputs,0,5);

	for(i = 0; i <= 4; i++)
	{
		modelOutput->SO[40*i][states[40*i]++] = i;
		modelOutput->OS[i][outputs[i]++] = 40*i;
	}
	simulator->model = QSS_Model(MOD_definition,MOD_dependencies,NULL,NULL,NULL);
	free(outputs);
	free(states);
}

void
CLC_initializeDataStructs (CLC_simulator simulator)
{
}
