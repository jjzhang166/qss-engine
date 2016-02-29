#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <common/utils.h>


#include <common/model.h>
#include <qss/qss_model.h>
#include <classic/classic_model.h>
static CLC_data modelData = NULL;

double __Ron = 0;
double __Roff = 0;
double __U = 0;
double __L = 0;
double __R = 0;
double __w = 0;

void
MOD_settings(SD_simulationSettings settings)
{
	 settings->debug = 0;
	 settings->parallel = FALSE;
	 settings->hybrid = TRUE;
	 settings->method = 7;
}

void
MOD_definition(double *x, double *d, double *alg, double t, double *dx)
{
	modelData->scalarEvaluations++;
		dx[0] = 1000.0*(__U*sin(__w*t)-x[0]);
		dx[1] = (x[0]-x[1]*(__R+d[(0)]))/__L;
}

void
MOD_zeroCrossing(int i, double *x, double *d, double *alg, double t, double *zc)
{
	modelData->zeroCrossings++;
	switch(i)
	{
		case 0:
			zc[0] = x[1]-(0.0);
			return;
		case 1:
			zc[0] = x[0]-(0.0);
			return;
	}
}

void
MOD_handlerPos(int i, double *x, double *d, double *alg, double t)
{
	switch(i)
	{
		case 1:
			d[(0)] = __Ron;
			return;
	}
}

void
MOD_handlerNeg(int i, double *x, double *d, double *alg, double t)
{
	switch(i)
	{
		case 0:
			d[(0)] = __Roff;
			return;
	}
}

void
MOD_output(int i, double *x, double *d, double *alg, double t, double *out)
{
	switch(i)
	{
		case 0:
			out[0] = x[1];
			return;
		case 1:
			out[0] = x[0];
			return;
	}
}

void
CLC_initializeDataStructs(CLC_simulator simulator)
{
	int discretes[1];
	int i = 0;
	int outputs[2];
	int states[2];
	simulator->data = CLC_Data(2,1,2,1,0,"rectifierDASSL");
modelData = simulator->data;

	// Allocate main data structures.
	__Ron = 1.000000000000000081803054e-05;
	__Roff = 1.000000000000000000000000e+05;
	__U = 311.0;
	__L = 1.000000000000000020816682e-03;
	__R = 10.0;
	__w = 3.141600000000000250111043e+02;
	modelData->d[(0)] = 1.000000000000000000000000e+05;
	// Initialize model code.
		modelData->event[0].direction = -1;
		modelData->event[1].direction = 1;
		double period[1];
	period[0] = 0.0002;
	simulator->output = SD_Output("rectifierDASSL",2,1,2,period,1,0,CI_Sampled,SD_Memory,MOD_output);
	SD_output modelOutput = simulator->output;

		modelOutput->nOS[0] = 1;
		modelOutput->nSO[1]++;
		modelOutput->nOS[1] = 1;
		modelOutput->nSO[0]++;
	SD_allocOutputMatrix(modelOutput,2,1);
		cleanVector(states,0,2);

		cleanVector(outputs,0,2);

		sprintf(modelOutput->variable[0].name,"iL");
		sprintf(modelOutput->variable[1].name,"u");
		cleanVector(outputs,0,2);

		modelOutput->SO[1][states[1]++] = 0;
		modelOutput->OS[0][outputs[0]++] = 1;
		modelOutput->SO[0][states[0]++] = 1;
		modelOutput->OS[1][outputs[1]++] = 0;
	simulator->model = CLC_Model(MOD_definition,MOD_zeroCrossing,MOD_handlerPos,MOD_handlerNeg);
}

void
QSS_initializeDataStructs (QSS_simulator simulator)
{
}
