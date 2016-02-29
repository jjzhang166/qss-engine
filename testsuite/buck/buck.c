#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <common/utils.h>


#include <common/model.h>
#include <qss/qss_model.h>
#include <classic/classic_model.h>

double __C;
double __L;
double __R;
double __U;
double __T;
double __DC;
double __ROn;
double __ROff;

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
	switch(i)
	{
		case 0:
			alg[0] = (x[0]*d[(1)]-__U)/(d[(1)]+d[(0)]);
			dx[1] = (-alg[0]*d[(0)]-x[3])/__L;
			return;
		case 1:
			dx[1] = (x[0]-x[3]/__R)/__C;
			return;
	}
}

void
MOD_dependencies(int i, double *x, double *d, double *alg, double t, double *der, int *map)
{
	switch(i)
	{
		case 0:
			alg[0] = (x[0]*d[(1)]-__U)/(d[(1)]+d[(0)]);
			der[0 + 1] = (-alg[0]*d[(0)]-x[3])/__L;
			der[3 + 1] = (x[0]-x[3]/__R)/__C;
			return;
		case 1:
			alg[0] = (x[0]*d[(1)]-__U)/(d[(1)]+d[(0)]);
			der[0 + 1] = (-alg[0]*d[(0)]-x[3])/__L;
			der[3 + 1] = (x[0]-x[3]/__R)/__C;
			return;
	}
}

void
MOD_zeroCrossing(int i, double *x, double *d, double *alg, double t, double *zc)
{
	switch(i)
	{
		case 0:
			zc[0] = t-(d[(2)]);
			return;
		case 1:
			zc[0] = t-d[(3)]-__DC*__T-(0.0);
			return;
		case 2:
			alg[0] = (x[0]*d[(1)]-__U)/(d[(1)]+d[(0)]);
			zc[0] = alg[0]*d[(4)]-(0.0);
			return;
	}
}

void
MOD_handlerPos(int i, double *x, double *d, double *alg, double t)
{
	switch(i)
	{
		case 0:
			d[(3)] = d[(2)];
			d[(2)] = d[(2)]+__T;
			d[(1)] = __ROn;
			return;
		case 1:
			d[(1)] = __ROff;
			d[(0)] = __ROn;
			d[(4)] = 1.0;
			return;
	}
}

void
MOD_handlerNeg(int i, double *x, double *d, double *alg, double t)
{
	switch(i)
	{
		case 2:
			d[(0)] = __ROff;
			d[(4)] = 0.0;
			return;
	}
}

void
MOD_output(int i, double *x, double *d, double *alg, double t, double *out)
{
	switch(i)
	{
		case 0:
			out[0] = x[3];
			return;
		case 1:
			out[0] = x[0];
			return;
	}
}

void
QSS_initializeDataStructs(QSS_simulator simulator)
{
	int *discretes = (int*)malloc(5*sizeof(int));
	int *events = (int*)malloc(3*sizeof(int));
	int *outputs = (int*)malloc(2*sizeof(int));
	int *states = (int*)malloc(2*sizeof(int));
	int i;
	simulator->data = QSS_Data(2,5,3,0,1,"buck");
QSS_data modelData = simulator->data;

	// Allocate main data structures.
	__C = 1.000000000000000047921736e-04;
	__L = 1.000000000000000047921736e-04;
	__R = 10.0;
	__U = 24.0;
	__T = 1.000000000000000047921736e-04;
	__DC = 5.000000000000000000000000e-01;
	__ROn = 1.000000000000000081803054e-05;
	__ROff = 1.000000000000000000000000e+05;
	// Initialize model code.
		modelData->d[(2)] = 0.0;
		modelData->d[(1)] = 1.000000000000000000000000e+05;
		modelData->d[(0)] = 1.000000000000000000000000e+05;
	modelData->nDS[0] = 1;
	modelData->nDS[1] = 2;
	modelData->nDS[0]++;
	modelData->nSD[1]++;
	modelData->nSD[0]++;
	modelData->nSD[1]++;
	modelData->nSD[0]++;
	modelData->nZS[2]++;
	modelData->nSZ[0]++;
	modelData->nHZ[0] = 2;
	modelData->nHZ[0] += 1;
	modelData->nHZ[1] = 1;
	modelData->nHZ[1] += 2;
	modelData->nHZ[2] = 1;
	modelData->nHZ[2] += 1;
	modelData->nHD[0] = 1;
	modelData->nHD[1] = 1;
	modelData->nHD[2] = 1;
	modelData->event[0].nLHSDsc = 3;
	modelData->event[1].nLHSDsc = 3;
	modelData->event[2].nLHSDsc = 2;
	QSS_allocDataMatrix(modelData);
	// Initialize model data.
	// Initialize model time.
		cleanVector(states,0,2);

	modelData->DS[0][states[0]++] = 1;
	modelData->DS[1][states[1]++] = 0;
	modelData->DS[1][states[1]++] = 1;
	modelData->DS[0][states[0]++] = 0;
		cleanVector(states,0,2);

	modelData->SD[1][states[1]++] = 0;
	modelData->SD[0][states[0]++] = 1;
	modelData->SD[1][states[1]++] = 1;
	modelData->SD[0][states[0]++] = 0;
		cleanVector(events,0,3);

	modelData->ZS[2][events[2]++] = 0;
		cleanVector(states,0,2);

	modelData->SZ[0][states[0]++] = 2;
		cleanVector(events,0,3);

	modelData->HZ[0][events[0]++] = 0;
	modelData->HZ[0][events[0]++] = 1;
	modelData->HZ[0][events[0]++] = 2;
	modelData->HZ[1][events[1]++] = 2;
	modelData->HZ[1][events[1]++] = 2;
	modelData->HZ[1][events[1]++] = 2;
	modelData->HZ[2][events[2]++] = 2;
	modelData->HZ[2][events[2]++] = 2;
		cleanVector(events,0,3);

	modelData->HD[1][events[1]++] = 0;
	modelData->HD[2][events[2]++] = 0;
	modelData->HD[0][events[0]++] = 0;
		cleanVector(events,0,3);

	modelData->event[0].LHSDsc[events[0]++] = 1;
	modelData->event[0].LHSDsc[events[0]++] = 2;
	modelData->event[0].LHSDsc[events[0]++] = 3;
	modelData->event[1].LHSDsc[events[1]++] = 0;
	modelData->event[1].LHSDsc[events[1]++] = 1;
	modelData->event[1].LHSDsc[events[1]++] = 4;
	modelData->event[2].LHSDsc[events[2]++] = 0;
	modelData->event[2].LHSDsc[events[2]++] = 4;
		cleanVector(events,0,3);

	modelData->event[0].direction = 1;
	modelData->event[0].relation = 2;
	modelData->event[1].direction = 1;
	modelData->event[1].relation = 2;
	modelData->event[2].direction = -1;
	modelData->event[2].relation = 0;
	simulator->time = QSS_Time(2,3,0,0,ST_Binary,NULL);

	simulator->output = SD_Output("buck",2,5,2,NULL,0,0,CI_Step,SD_Memory,MOD_output);
SD_output modelOutput = simulator->output;

		modelOutput->nOS[0] = 1;
		modelOutput->nSO[1]++;
		modelOutput->nOS[1] = 1;
		modelOutput->nSO[0]++;
	SD_allocOutputMatrix(modelOutput,2,5);
		cleanVector(states,0,2);

		cleanVector(outputs,0,2);

		sprintf(modelOutput->variable[0].name,"uC");
		sprintf(modelOutput->variable[1].name,"iL");
		cleanVector(outputs,0,2);

		modelOutput->SO[1][states[1]++] = 0;
		modelOutput->OS[0][outputs[0]++] = 1;
		modelOutput->SO[0][states[0]++] = 1;
		modelOutput->OS[1][outputs[1]++] = 0;
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
