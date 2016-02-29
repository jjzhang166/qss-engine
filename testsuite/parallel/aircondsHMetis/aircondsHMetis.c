#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "pkg_math.h"
#include <common/utils.h>


#include <common/model.h>
#include <qss/qss_model.h>
#include <classic/classic_model.h>

double __CAP[5000];
double __RES[5000];
double __POT[5000];
double __THA = 0;
double __pmax = 0;

void
MOD_settings(SD_simulationSettings settings)
{
	 settings->debug = 0;
	 settings->parallel = TRUE;
	 settings->hybrid = TRUE;
	 settings->method = 5;
}

void
MOD_definition(int i, double *x, double *d, double *alg, double t, double *dx)
{
	int j = 0;
	j = i;
	if(j >=0 && j <= 4999)
	{
		dx[1] = (__THA/__RES[(j)]-__POT[(j)]*d[(j+1)]-x[(j) * 4]/__RES[(j)]+d[(j+15001)]/__RES[(j)])/__CAP[(j)];
		dx[2] = (-(1.0/(__RES[(j)]))*x[(j) * 4 + 1]*(1.0/(__CAP[(j)])))/2;
		dx[3] = (-(1.0/(__RES[(j)]))*(1.0/(__CAP[(j)]))*x[(j) * 4 + 2]*2)/6;
	}
}

void
MOD_dependencies(int i, double *x, double *d, double *alg, double t, double *der, int *map)
{
	int j = 0;
	j = i;
	if(j >=0 && j <= 4999)
	{
	if (map[j] != NOT_ASSIGNED)
{
		der[(j) * 4 + 1] = (__THA/__RES[(j)]-__POT[(j)]*d[(j+1)]-x[(j) * 4]/__RES[(j)]+d[(j+15001)]/__RES[(j)])/__CAP[(j)];
		der[(j) * 4 + 2] = (-(1.0/(__RES[(j)]))*x[(j) * 4 + 1]*(1.0/(__CAP[(j)])))/2;
		der[(j) * 4 + 3] = (-(1.0/(__RES[(j)]))*(1.0/(__CAP[(j)]))*x[(j) * 4 + 2]*2)/6;
	}
	}
}

void
MOD_zeroCrossing(int i, double *x, double *d, double *alg, double t, double *zc)
{
	if(i >= 0 && i <= 4999)
	{
		zc[0] = x[(i) * 4]-d[(i+5001)]+d[(i+1)]-5.000000000000000000000000e-01-(0.0);
		zc[1] = x[(i) * 4 + 1];
		zc[2] = (x[(i) * 4 + 2]*2)/2;
	}
	if(i >= 5000 && i <= 9999)
	{
		zc[0] = t-(1000.0);
		zc[1] = 1.0;
		zc[2] = (0.0)/2;
	}
	if(i >= 10000 && i <= 14999)
	{
		zc[0] = t-(2000.0);
		zc[1] = 1.0;
		zc[2] = (0.0)/2;
	}
	if(i >= 15000 && i <= 19999)
	{
		zc[0] = t-(d[(i-4999)]);
		zc[1] = 1.0;
		zc[2] = (0.0)/2;
	}
}

void
MOD_handlerPos(int i, double *x, double *d, double *alg, double t)
{
	if(i >= 0 && i <= 4999)
	{
		d[(i+1)] = 1.0;
		d[(0)] = d[(0)]+__POT[(i)];
	}
	if(i >= 5000 && i <= 9999)
	{
		d[(i+1)] = 2.050000000000000000000000e+01;
	}
	if(i >= 10000 && i <= 14999)
	{
		d[(i-4999)] = 20.0;
	}
	if(i >= 15000 && i <= 19999)
	{
		d[(i-4999)] = d[(i-4999)]+1.0;
		d[(i+1)] = __math__rand(2.0)-1.0;
	}
}

void
MOD_handlerNeg(int i, double *x, double *d, double *alg, double t)
{
	if(i >= 0 && i <= 4999)
	{
	if(t>0.0)
	{
		d[(i+1)] = 0.0;
		d[(0)] = d[(0)]-__POT[(i)];
	}
	}
}

void
MOD_output(int i, double *x, double *d, double *alg, double t, double *out)
{
	switch(i)
	{
		case 0:
			out[0] = d[(0)];
			return;
	}
}

void
QSS_initializeDataStructs(QSS_simulator simulator)
{
	int *discretes = (int*)malloc(20001*sizeof(int));
	int *events = (int*)malloc(20000*sizeof(int));
	int *outputs = (int*)malloc(1*sizeof(int));
	int *states = (int*)malloc(5000*sizeof(int));
	int i2;
	int i;
	int j = 0;
	simulator->data = QSS_Data(5000,20001,20000,0,0,"aircondsHMetis");
QSS_data modelData = simulator->data;

	// Allocate main data structures.
	__THA = 32.0;
	__pmax = 0.0;
	// Initialize model code.
	for(i2 = 0; i2 <= 4999; i2++)
	{
		modelData->x[(i2) * 4] = __math__rand(4.0)+18.0;
		__CAP[(i2)] = __math__rand(100.0)+550.0;
		__RES[(i2)] = __math__rand(4.000000000000000222044605e-01)+1.800000000000000044408921e+00;
		__POT[(i2)] = __math__rand(2.0)+13.0;
		__pmax = __pmax+__POT[(i2)];
		modelData->d[(i2+10001)] = 1.0;
		modelData->d[(i2+15001)] = __math__rand(2.0)-1.0;
		modelData->d[(i2+5001)] = 20.0;
	}
	for(i = 0; i <= 4999; i++)
	{
		modelData->nDS[i] = 1;
	}
	for(i = 0; i <= 4999; i++)
	{
		modelData->nSD[i]++;
	}
	for(i = 0; i <= 4999; i++)
	{
		modelData->nZS[i] = 1;
	}
	for(i = 0; i <= 4999; i++)
	{
		modelData->nSZ[i]++;
	}
	for(i = 0; i <= 4999; i++)
	{
		modelData->nHZ[i]++;
	}
	for(i = 5000; i <= 9999; i++)
	{
		modelData->nHZ[i]++;
	}
	for(i = 10000; i <= 14999; i++)
	{
		modelData->nHZ[i]++;
	}
	for(i = 15000; i <= 19999; i++)
	{
		modelData->nHZ[i]++;
	}
	for(i = 0; i <= 4999; i++)
	{
		modelData->nHD[i]++;
	}
	for(i = 15000; i <= 19999; i++)
	{
		modelData->nHD[i]++;
	}
	for(i = 0; i <= 4999; i++)
	{
		modelData->event[i].nLHSDsc = 2;
		modelData->event[i+5000].nLHSDsc = 1;
		modelData->event[i+10000].nLHSDsc = 1;
		modelData->event[i+15000].nLHSDsc = 2;
	}
	QSS_allocDataMatrix(modelData);
	// Initialize model data.
	// Initialize model time.
		cleanVector(states,0,5000);

	for(i = 0; i <= 4999; i++)
	{
		modelData->DS[i][states[i]++] = i;
	}
		cleanVector(states,0,5000);

	for(i = 0; i <= 4999; i++)
	{
		modelData->SD[i][states[i]++] = i;
	}
		cleanVector(events,0,20000);

	for(i = 0; i <= 4999; i++)
	{
		modelData->ZS[i][events[i]++] = i;
	}
		cleanVector(states,0,5000);

	for(i = 0; i <= 4999; i++)
	{
		modelData->SZ[i][states[i]++] = i;
	}
		cleanVector(events,0,20000);

	for(i = 0; i <= 4999; i++)
	{
		modelData->HZ[i][events[i]++] = i;
	}
	for(i = 5000; i <= 9999; i++)
	{
		modelData->HZ[i][events[i]++] = i-5000;
	}
	for(i = 10000; i <= 14999; i++)
	{
		modelData->HZ[i][events[i]++] = i-10000;
	}
	for(i = 15000; i <= 19999; i++)
	{
		modelData->HZ[i][events[i]++] = i;
	}
		cleanVector(events,0,20000);

	for(i = 0; i <= 4999; i++)
	{
		modelData->HD[i][events[i]++] = i;
	}
	for(i = 15000; i <= 19999; i++)
	{
		modelData->HD[i][events[i]++] = i-15000;
	}
		cleanVector(events,0,20000);

	for(i = 0; i <= 4999; i++)
	{
		modelData->event[i].LHSDsc[events[i]++] = 0;
		modelData->event[i].LHSDsc[events[i]++] = i+1;
		modelData->event[i+5000].LHSDsc[events[i+5000]++] = i+5001;
		modelData->event[i+10000].LHSDsc[events[i+10000]++] = i+5001;
		modelData->event[i+15000].LHSDsc[events[i+15000]++] = i+10001;
		modelData->event[i+15000].LHSDsc[events[i+15000]++] = i+15001;
	}
		cleanVector(events,0,20000);

	for(i = 0; i <= 4999; i++)
	{
		modelData->event[i].direction = 0;
		modelData->event[i].relation = 2;
		modelData->event[i+5000].direction = 1;
		modelData->event[i+5000].relation = 2;
		modelData->event[i+10000].direction = 1;
		modelData->event[i+10000].relation = 2;
		modelData->event[i+15000].direction = 1;
		modelData->event[i+15000].relation = 2;
	}
	simulator->time = QSS_Time(5000,20000,0,0,ST_Binary,NULL);

	simulator->output = SD_Output("aircondsHMetis",1,20001,5000,NULL,0,0,CI_Step,SD_Memory,MOD_output);
SD_output modelOutput = simulator->output;

		modelOutput->nOD[0] = 1;
		modelOutput->nDO[0]++;
	SD_allocOutputMatrix(modelOutput,5000,20001);
		cleanVector(discretes,0,20001);

		cleanVector(outputs,0,1);

		sprintf(modelOutput->variable[0].name,"ptotal");
		cleanVector(outputs,0,1);

		modelOutput->DO[0][discretes[0]++] = 0;
		modelOutput->OD[0][outputs[0]++] = 0;
	simulator->model = QSS_Model(MOD_definition,MOD_dependencies,MOD_zeroCrossing,MOD_handlerPos,MOD_handlerNeg);
	free(discretes);
	free(events);
	free(outputs);
	free(states);
}

void
CLC_initializeDataStructs (CLC_simulator simulator)
{
}
