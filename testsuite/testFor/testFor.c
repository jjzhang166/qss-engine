#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <common/utils.h>


#include <common/model.h>
#include <qss/qss_model.h>
#include <classic/classic_model.h>

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
			dx[1] = (-x[0]*__ax/__dx)+(-x[0]*__ay/__dy);
			return;
		default:
			j = i;
			if(j >=1 && j <= 1)
			{
				dx[1] = (-x[(j) * 3]+x[(j-1) * 3])*__ax/__dx+(-x[(j) * 3]*__ay/__dy);
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
			der[0 + 1] = (-x[0]*__ax/__dx)+(-x[0]*__ay/__dy);
			break;
	}
	j = i+1;
	if(j >=1 && j <= 1)
	{
		der[(j) * 3 + 1] = (-x[(j) * 3]+x[(j-1) * 3])*__ax/__dx+(-x[(j) * 3]*__ay/__dy);
	}
	j = i;
	if(j >=1 && j <= 1)
	{
		der[(j) * 3 + 1] = (-x[(j) * 3]+x[(j-1) * 3])*__ax/__dx+(-x[(j) * 3]*__ay/__dy);
	}
}

void
MOD_output(int i, double *x, double *d, double *alg, double t, double *out)
{
	int j = 0;
	j = i;
	if(j >=0 && j <= 1)
	{
		out[0] = x[(j) * 3];
	}
}

void
QSS_initializeDataStructs(QSS_simulator simulator)
{
	int *outputs = (int*)malloc(2*sizeof(int));
	int *states = (int*)malloc(2*sizeof(int));
	int i0;
	int i;
	int j = 0;
	simulator->data = QSS_Data(2,0,0,0,0,"testFor");
QSS_data modelData = simulator->data;

	// Allocate main data structures.
	__ax = 1.0;
	__ay = 1.0;
	__r = 1000.0;
	__dx = 10.0/100;
	__dy = 10.0/100;
	// Initialize model code.
	for(i0 = 0; i0 <= 1; i0++)
	{
		modelData->x[(i0) * 3] = 1.0;
	}
	modelData->nDS[0] = 1;
	for(i = 1; i <= 1; i++)
	{
		modelData->nDS[i] = 2;
	}
	modelData->nSD[0]++;
	for(i = 1; i <= 1; i++)
	{
		modelData->nSD[i-1]++;
		modelData->nSD[i]++;
	}
	QSS_allocDataMatrix(modelData);
	// Initialize model data.
	// Initialize model time.
		cleanVector(states,0,2);

	modelData->DS[0][states[0]++] = 0;
	for(i = 1; i <= 1; i++)
	{
		modelData->DS[i][states[i]++] = i-1;
		modelData->DS[i][states[i]++] = i;
	}
		cleanVector(states,0,2);

	modelData->SD[0][states[0]++] = 0;
	for(i = 1; i <= 1; i++)
	{
		modelData->SD[i-1][states[i-1]++] = i;
		modelData->SD[i][states[i]++] = i;
	}
	simulator->time = QSS_Time(2,0,0,0,ST_Binary,NULL);

	simulator->output = SD_Output("testFor",2,0,2,NULL,0,0,CI_Step,SD_Memory,MOD_output);
SD_output modelOutput = simulator->output;

	for(i = 0; i <= 1; i++)
	{
		modelOutput->nOS[i] = 1;
		modelOutput->nSO[i]++;
	}
	SD_allocOutputMatrix(modelOutput,2,0);
		cleanVector(states,0,2);

		cleanVector(outputs,0,2);

	for(i = 0; i <= 1; i++)
	{
		sprintf(modelOutput->variable[i].name,"u1[%d]",i+1);
	}
		cleanVector(outputs,0,2);

	for(i = 0; i <= 1; i++)
	{
		modelOutput->SO[i][states[i]++] = i;
		modelOutput->OS[i][outputs[i]++] = i;
	}
	simulator->model = QSS_Model(MOD_definition,MOD_dependencies,NULL,NULL,NULL);
	free(outputs);
	free(states);
}

void
CLC_initializeDataStructs (CLC_simulator simulator)
{
}
