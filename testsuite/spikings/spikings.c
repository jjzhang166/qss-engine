#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "pkg_math.h"
#include <common/utils.h>


#include <common/model.h>
#include <qss/qss_model.h>
#include <classic/classic_model.h>

double __a;
double __b;
double __c[1000];
double __d[1000];
double __I[1000];
double __r[1000];

void
MOD_settings(SD_simulationSettings settings)
{
	 settings->debug = 0;
	 settings->parallel = FALSE;
	 settings->hybrid = TRUE;
	 settings->method = 4;
}

void
MOD_definition(int i, double *x, double *d, double *alg, double t, double *dx)
{
	int j = 0;
	j = i;
	if(j >=0 && j <= 999)
	{
		dx[1] = 4.000000000000000083266727e-02*x[(j) * 3]*x[(j) * 3]+5.0*x[(j) * 3]+140.0-x[(j+1000) * 3]+__I[(j)];
	}
	j = i-1000;
	if(j >=0 && j <= 999)
	{
		dx[1] = __a*(__b*x[(j) * 3]-x[(j+1000) * 3]);
	}
}

void
MOD_dependencies(int i, double *x, double *d, double *alg, double t, double *der, int *map)
{
	int j = 0;
	j = i;
	if(j >=0 && j <= 999)
	{
		der[(j) * 3 + 1] = 4.000000000000000083266727e-02*x[(j) * 3]*x[(j) * 3]+5.0*x[(j) * 3]+140.0-x[(j+1000) * 3]+__I[(j)];
		der[(j+1000) * 3 + 1] = __a*(__b*x[(j) * 3]-x[(j+1000) * 3]);
	}
	j = i-1000;
	if(j >=0 && j <= 999)
	{
		der[(j) * 3 + 1] = 4.000000000000000083266727e-02*x[(j) * 3]*x[(j) * 3]+5.0*x[(j) * 3]+140.0-x[(j+1000) * 3]+__I[(j)];
		der[(j+1000) * 3 + 1] = __a*(__b*x[(j) * 3]-x[(j+1000) * 3]);
	}
}

void
MOD_zeroCrossing(int i, double *x, double *d, double *alg, double t, double *zc)
{
	if(i >= 0 && i <= 999)
	{
		zc[0] = x[(i) * 3]-(30.0);
	}
}

void
MOD_handlerPos(int i, double *x, double *d, double *alg, double t)
{
	if(i >= 0 && i <= 999)
	{
		x[(i) * 3] = __c[(i)];
		x[(i+1000) * 3] = x[(i+1000) * 3]+__d[(i)];
	}
}

void
MOD_output(int i, double *x, double *d, double *alg, double t, double *out)
{
	switch(i)
	{
		case 0:
			out[0] = x[0];
			return;
		case 1:
			out[0] = x[3];
			return;
	}
}

void
QSS_initializeDataStructs(QSS_simulator simulator)
{
	int *events = (int*)malloc(1000*sizeof(int));
	int *outputs = (int*)malloc(2*sizeof(int));
	int *states = (int*)malloc(2000*sizeof(int));
	int i0;
	int i;
	int j = 0;
	simulator->data = QSS_Data(2000,0,1000,0,0,"spikings");
QSS_data modelData = simulator->data;

	// Allocate main data structures.
	__a = 2.000000000000000041633363e-02;
	__b = 2.000000000000000111022302e-01;
	// Initialize model code.
	for(i0 = 0; i0 <= 999; i0++)
	{
		modelData->x[(i0) * 3] = -60.0;
		modelData->x[(i0+1000) * 3] = 0.0;
		__r[(i0)] = __math__rand(1.0);
		__c[(i0)] = -65.0+15.0*pow(__r[(i0)],2.0);
		__d[(i0)] = 8.0-6.0*pow(__r[(i0)],2.0);
		__I[(i0)] = __math__rand(10.0);
	}
	for(i = 0; i <= 999; i++)
	{
		modelData->nDS[i] = 2;
		modelData->nDS[i+1000] = 2;
	}
	for(i = 0; i <= 999; i++)
	{
		modelData->nSD[i]++;
		modelData->nSD[i+1000]++;
		modelData->nSD[i]++;
		modelData->nSD[i+1000]++;
	}
	for(i = 0; i <= 999; i++)
	{
		modelData->nZS[i] = 1;
	}
	for(i = 0; i <= 999; i++)
	{
		modelData->nSZ[i]++;
	}
	for(i = 0; i <= 999; i++)
	{
		modelData->nHZ[i]++;
	}
	for(i = 0; i <= 999; i++)
	{
		modelData->event[i].nLHSSt = 2;
	}
	for(i = 0; i <= 999; i++)
	{
		modelData->event[i].nRHSSt = 1;
	}
	QSS_allocDataMatrix(modelData);
	// Initialize model data.
	// Initialize model time.
		cleanVector(states,0,2000);

	for(i = 0; i <= 999; i++)
	{
		modelData->DS[i][states[i]++] = i;
		modelData->DS[i][states[i]++] = i+1000;
		modelData->DS[i+1000][states[i+1000]++] = i;
		modelData->DS[i+1000][states[i+1000]++] = i+1000;
	}
		cleanVector(states,0,2000);

	for(i = 0; i <= 999; i++)
	{
		modelData->SD[i][states[i]++] = i;
		modelData->SD[i+1000][states[i+1000]++] = i;
		modelData->SD[i][states[i]++] = i+1000;
		modelData->SD[i+1000][states[i+1000]++] = i+1000;
	}
		cleanVector(events,0,1000);

	for(i = 0; i <= 999; i++)
	{
		modelData->ZS[i][events[i]++] = i;
	}
		cleanVector(states,0,2000);

	for(i = 0; i <= 999; i++)
	{
		modelData->SZ[i][states[i]++] = i;
	}
		cleanVector(events,0,1000);

	for(i = 0; i <= 999; i++)
	{
		modelData->HZ[i][events[i]++] = i;
	}
		cleanVector(events,0,1000);

	for(i = 0; i <= 999; i++)
	{
		modelData->event[i].LHSSt[events[i]++] = i;
		modelData->event[i].LHSSt[events[i]++] = i+1000;
	}
		cleanVector(events,0,1000);

	for(i = 0; i <= 999; i++)
	{
		modelData->event[i].RHSSt[events[i]++] = i+1000;
	}
		cleanVector(events,0,1000);

	for(i = 0; i <= 999; i++)
	{
		modelData->event[i].direction = 1;
		modelData->event[i].relation = 2;
	}
	simulator->time = QSS_Time(2000,1000,0,0,ST_Binary,NULL);

	simulator->output = SD_Output("spikings",2,0,2000,NULL,0,0,CI_Step,SD_Memory,MOD_output);
SD_output modelOutput = simulator->output;

		modelOutput->nOS[0] = 1;
		modelOutput->nSO[0]++;
		modelOutput->nOS[1] = 1;
		modelOutput->nSO[1]++;
	SD_allocOutputMatrix(modelOutput,2000,0);
		cleanVector(states,0,2000);

		cleanVector(outputs,0,2);

		sprintf(modelOutput->variable[0].name,"v[1]");
		sprintf(modelOutput->variable[1].name,"v[2]");
		cleanVector(outputs,0,2);

		modelOutput->SO[0][states[0]++] = 0;
		modelOutput->OS[0][outputs[0]++] = 0;
		modelOutput->SO[1][states[1]++] = 1;
		modelOutput->OS[1][outputs[1]++] = 1;
	simulator->model = QSS_Model(MOD_definition,MOD_dependencies,MOD_zeroCrossing,MOD_handlerPos,NULL);
	free(events);
	free(outputs);
	free(states);
}

void
CLC_initializeDataStructs (CLC_simulator simulator)
{
}
