#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "pkg_math.h"
#include <common/utils.h>


#include <common/model.h>
#include <common/commands.h>
#include <qss/qss_model.h>
#include <classic/classic_model.h>

double __PAR_CAP[200];
double __PAR_RES[200];
double __PAR_POT[200];
double __PAR_THA = 0;
double __PAR_pmax = 0;

void
MOD_settings(SD_simulationSettings settings)
{
	 settings->debug = 0;
	 settings->parallel = FALSE;
	 settings->hybrid = TRUE;
	 settings->method = 5;
}

void
MOD_definition(int i, double *x, double *d, double *alg, double t, double *dx)
{
	int j = 0;
	j = i;
	if(j >=0 && j <= 199)
	{
		dx[1] = (__PAR_THA/__PAR_RES[(j)]-__PAR_POT[(j)]*d[(j+1)]-x[(j) * 4]/__PAR_RES[(j)]+d[(j+601)]/__PAR_RES[(j)])/__PAR_CAP[(j)];
	}
}

void
MOD_dependencies(int i, double *x, double *d, double *alg, double t, double *der, int *map)
{
	int j = 0;
	j = i;
	if(j >=0 && j <= 199)
	{
		der[(j) * 4 + 1] = (__PAR_THA/__PAR_RES[(j)]-__PAR_POT[(j)]*d[(j+1)]-x[(j) * 4]/__PAR_RES[(j)]+d[(j+601)]/__PAR_RES[(j)])/__PAR_CAP[(j)];
	}
}

void
MOD_zeroCrossing(int i, double *x, double *d, double *alg, double t, double *zc)
{
	if(i >= 0 && i <= 199)
	{
		zc[0] = x[(i) * 4]-d[(i+201)]+d[(i+1)]-5.000000000000000000000000e-01-(0.0);
	}
	if(i >= 200 && i <= 399)
	{
		zc[0] = t-(1000.0);
	}
	if(i >= 400 && i <= 599)
	{
		zc[0] = t-(2000.0);
	}
	if(i >= 600 && i <= 799)
	{
		zc[0] = t-(d[(i-199)]);
	}
}

void
MOD_handlerPos(int i, double *x, double *d, double *alg, double t)
{
	if(i >= 0 && i <= 199)
	{
		d[(i+1)] = 1.0;
		d[(0)] = d[(0)]+__PAR_POT[(i)];
	}
	if(i >= 200 && i <= 399)
	{
		d[(i+1)] = 2.050000000000000000000000e+01;
	}
	if(i >= 400 && i <= 599)
	{
		d[(i-199)] = 20.0;
	}
	if(i >= 600 && i <= 799)
	{
		d[(i-199)] = d[(i-199)]+1.0;
		d[(i+1)] = __math__rand(2.0)-1.0;
	}
}

void
MOD_handlerNeg(int i, double *x, double *d, double *alg, double t)
{
	if(i >= 0 && i <= 199)
	{
		d[(i+1)] = 0.0;
		d[(0)] = d[(0)]-__PAR_POT[(i)];
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
	int *discretes = (int*)malloc(801*sizeof(int));
	int *events = (int*)malloc(800*sizeof(int));
	int *outputs = (int*)malloc(1*sizeof(int));
	int *states = (int*)malloc(200*sizeof(int));
	int i2;
	int i;
	int j = 0;
	simulator->data = QSS_Data(200,801,800,0,0,"airconds");
QSS_data modelData = simulator->data;

	// Allocate main data structures.
	__PAR_THA = 32.0;
	__PAR_pmax = 0.0;
	// Initialize model code.
	for(i2 = 0; i2 <= 199; i2++)
	{
		modelData->x[(i2) * 4] = __math__rand(4.0)+18.0;
		__PAR_CAP[(i2)] = __math__rand(100.0)+550.0;
		__PAR_RES[(i2)] = __math__rand(4.000000000000000222044605e-01)+1.800000000000000044408921e+00;
		__PAR_POT[(i2)] = __math__rand(2.0)+13.0;
		__PAR_pmax = __PAR_pmax+__PAR_POT[(i2)];
		modelData->d[(i2+401)] = 1.0;
		modelData->d[(i2+601)] = __math__rand(2.0)-1.0;
		modelData->d[(i2+201)] = 20.0;
	if(modelData->x[(i2) * 4]-modelData->d[(i2+201)]+modelData->d[(i2+1)]-5.000000000000000000000000e-01>0.0)
	{
		modelData->d[(i2+1)] = 1.0;
		modelData->d[(0)] = modelData->d[(0)]+__PAR_POT[(i2)];
	}
	}
	for(i = 0; i <= 199; i++)
	{
		modelData->nDS[i] = 1;
	}
	for(i = 0; i <= 199; i++)
	{
		modelData->nSD[i]++;
	}
	for(i = 0; i <= 199; i++)
	{
		modelData->nZS[i] = 1;
	}
	for(i = 0; i <= 199; i++)
	{
		modelData->nSZ[i]++;
	}
	for(i = 0; i <= 199; i++)
	{
		modelData->nHZ[i]++;
	}
	for(i = 200; i <= 399; i++)
	{
		modelData->nHZ[i]++;
	}
	for(i = 400; i <= 599; i++)
	{
		modelData->nHZ[i]++;
	}
	for(i = 600; i <= 799; i++)
	{
		modelData->nHZ[i]++;
	}
	for(i = 0; i <= 199; i++)
	{
		modelData->nHD[i]++;
	}
	for(i = 600; i <= 799; i++)
	{
		modelData->nHD[i]++;
	}
	for(i = 0; i <= 199; i++)
	{
		modelData->event[i].nLHSDsc = 2;
		modelData->event[i+200].nLHSDsc = 1;
		modelData->event[i+400].nLHSDsc = 1;
		modelData->event[i+600].nLHSDsc = 2;
	}
	QSS_allocDataMatrix(modelData);
	// Initialize model data.
	// Initialize model time.
		cleanVector(states,0,200);

	for(i = 0; i <= 199; i++)
	{
		modelData->DS[i][states[i]++] = i;
	}
		cleanVector(states,0,200);

	for(i = 0; i <= 199; i++)
	{
		modelData->SD[i][states[i]++] = i;
	}
		cleanVector(events,0,800);

	for(i = 0; i <= 199; i++)
	{
		modelData->ZS[i][events[i]++] = i;
	}
		cleanVector(states,0,200);

	for(i = 0; i <= 199; i++)
	{
		modelData->SZ[i][states[i]++] = i;
	}
		cleanVector(events,0,800);

	for(i = 0; i <= 199; i++)
	{
		modelData->HZ[i][events[i]++] = i;
	}
	for(i = 200; i <= 399; i++)
	{
		modelData->HZ[i][events[i]++] = i-200;
	}
	for(i = 400; i <= 599; i++)
	{
		modelData->HZ[i][events[i]++] = i-400;
	}
	for(i = 600; i <= 799; i++)
	{
		modelData->HZ[i][events[i]++] = i;
	}
		cleanVector(events,0,800);

	for(i = 0; i <= 199; i++)
	{
		modelData->HD[i][events[i]++] = i;
	}
	for(i = 600; i <= 799; i++)
	{
		modelData->HD[i][events[i]++] = i-600;
	}
		cleanVector(events,0,800);

	for(i = 0; i <= 199; i++)
	{
		modelData->event[i].LHSDsc[events[i]++] = 0;
		modelData->event[i].LHSDsc[events[i]++] = i+1;
		modelData->event[i+200].LHSDsc[events[i+200]++] = i+201;
		modelData->event[i+400].LHSDsc[events[i+400]++] = i+201;
		modelData->event[i+600].LHSDsc[events[i+600]++] = i+401;
		modelData->event[i+600].LHSDsc[events[i+600]++] = i+601;
	}
		cleanVector(events,0,800);

	for(i = 0; i <= 199; i++)
	{
		modelData->event[i].direction = 0;
		modelData->event[i].relation = 2;
		modelData->event[i+200].direction = 1;
		modelData->event[i+200].relation = 2;
		modelData->event[i+400].direction = 1;
		modelData->event[i+400].relation = 2;
		modelData->event[i+600].direction = 1;
		modelData->event[i+600].relation = 2;
	}
	simulator->time = QSS_Time(200,800,0,0,ST_Binary,NULL);

	simulator->output = SD_Output("airconds",1,801,200,NULL,0,0,CI_Step,SD_Memory,MOD_output);
SD_output modelOutput = simulator->output;

		modelOutput->nOD[0] = 1;
		modelOutput->nDO[0]++;
	SD_allocOutputMatrix(modelOutput,200,801);
		cleanVector(discretes,0,801);

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
