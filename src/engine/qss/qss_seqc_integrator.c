/*****************************************************************************

    This file is part of QSS Solver.

    QSS Solver is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    QSS Solver is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with QSS Solver.  If not, see <http://www.gnu.org/licenses/>.

******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>

#include <qss/qss_integrator.h>
#include <qss/qss_seqc_integrator.h>
#include <qss/qss_simulator.h>
#include <qss/qss_model.h>

void
QSS_SEQC_integrate (SIM_simulator simulate)
{
  QSS_simulator simulator = simulate->state->sim;
  int i, j;
  double elapsed;
  // Local data structure mappings.
  QSS_data qssData = simulator->data;
  QSS_time qssTime = simulator->time;
  FRW_framework frw = simulator->frw;
  OUT_output log = simulator->log;
  SC_scheduler scheduler = simulator->scheduler;
  QSS_model qssModel = simulator->model;
  QA_quantizer quantizer = simulator->quantizer;
  SD_output output = simulator->output;
#ifdef DEBUG
  SD_simulationSettings settings = simulator->settings;
  SD_simulationLog simulationLog = simulator->simulationLog;
#endif
  double t = qssTime->time;
  int index = qssTime->minIndex;
  QSS_StepType type = qssTime->type;
  int nSD, nOutputs = output->outputs;
  const double ft = qssData->ft;
  const int xOrder = qssData->order;
  const int coeffs = xOrder + 1;
  double *tq = qssTime->tq;
  double *tx = qssTime->tx;
  double *nextStateTime = qssTime->nextStateTime;
  double *dQRel = qssData->dQRel;
  double *dQMin = qssData->dQMin;
  double *lqu = qssData->lqu;
  double *x = qssData->x;
  double *q = qssData->q;
  int **SD = qssData->SD;
  int *TD = qssData->TD;
  int cf0, infCf0;
  getTime (simulator->iTime);
#ifdef SYNC_RT
  setInitRealTime();
#endif
#ifdef DEBUG
  if (settings->debug & SD_DBG_StepInfo)
    {
      SD_print (simulator->simulationLog, "Begin Simulation:");
    }
#endif
  while (t < ft)
    {
#ifdef SYNC_RT
      /* Sync */
      waitUntil(t);
#endif
#ifdef DEBUG
      if (settings->debug & SD_DBG_StepInfo)
	{
	  SD_print (simulationLog, "Simulation Time: %g", t);
	}
#endif
      switch (type)
	{
	case ST_State:
	  {
#ifdef DEBUG
	    if (settings->debug & SD_DBG_StepInfo)
	      {
		SD_print (simulationLog, "State Variable: %d", index);
	      }
	    if (settings->debug & SD_DBG_VarChanges)
	      {
		simulationLog->states[index]++;
	      }
#endif
	    cf0 = index * coeffs;
	    elapsed = t - tx[index];
	    advanceTime (cf0, elapsed, x, xOrder);
	    tx[index] = t;
	    lqu[index] = dQRel[index] * fabs (x[cf0]);
	    if (lqu[index] < dQMin[index])
	      {
		lqu[index] = dQMin[index];
	      }
	    QA_updateQuantizedState (quantizer, index, q, x, lqu);
	    tq[index] = t;
	    QA_nextTime (quantizer, index, t, nextStateTime, x, lqu);
	    nSD = qssData->nSD[index];
	    for (i = 0; i < nSD; i++)
	      {
		j = SD[index][i];
		elapsed = t - tx[j];
		infCf0 = j * coeffs;
		if (elapsed > 0)
		  {
		    x[infCf0] = evaluatePoly (infCf0, elapsed, x, xOrder);
		    tx[j] = t;
		  }
	      }
	    FRW_recomputeDerivatives (frw, qssModel, qssData, qssTime, index);
	    QA_recomputeNextTimes (quantizer, nSD, SD[index], t, nextStateTime,
				   x, lqu, q);
	    if (nOutputs)
	      {
		if (output->nSO[index])
		  {
		    OUT_write (log, qssData, qssTime, output);
		  }
	      }
	  }
	  break;
	case ST_Input:
	  {
#ifdef DEBUG
	    if (settings->debug & SD_DBG_StepInfo)
	      {
		SD_print (simulationLog, "Input: %d", index);
	      }
#endif
	    j = TD[index];
	    elapsed = t - tx[j];
	    infCf0 = j * coeffs;
	    if (elapsed > 0)
	      {
		x[infCf0] = evaluatePoly (infCf0, elapsed, x, xOrder);
		tx[j] = t;
	      }
	    FRW_recomputeDerivative (frw, qssModel, qssData, qssTime, j);
	    FRW_nextInputTime (frw, qssModel, qssData, qssTime, elapsed, j,
			       index);
	    QA_recomputeNextTime (quantizer, j, t, nextStateTime, x, lqu, q);
	  }
	  break;
	default:
	  break;
	}
      SC_update (scheduler, qssData, qssTime);
      t = qssTime->time;
      index = qssTime->minIndex;
      type = qssTime->type;
      simulator->totalSteps++;
#ifdef DEBUG
      if (settings->debug & SD_DBG_StepInfo)
	{
	  SD_print (simulationLog, "");
	}
#endif
    }
  getTime (simulator->sTime);
  subTime (simulator->sTime, simulator->iTime);
  simulator->simulationTime = getTimeValue (simulator->sTime);
  QSS_SEQ_saveLog (simulator);
  QSS_SEQ_printSimulationLog (simulator);
}
