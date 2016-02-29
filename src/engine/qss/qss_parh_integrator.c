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

#include <qss/qss_parh_integrator.h>

#ifdef __linux__
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>

#include <common/data.h>
#include <qss/qss_data.h>
#include <qss/qss_integrator.h>
#include <qss/qss_simulator.h>
#include <qss/qss_model.h>
#include <qss/qss_partition.h>
#include <qss/qss_parallel.h>
#include <qss/qss_lp.h>

void
QSS_PARH_externalEvent (QSS_simulator simulator, IBX_message message)
{
  int i, j;
  double elapsed;
  QSS_data qssData = simulator->data;
  QSS_time qssTime = simulator->time;
  FRW_framework frw = simulator->frw;
  QSS_model qssModel = simulator->model;
  QA_quantizer quantizer = simulator->quantizer;
  QSS_LP_data lp = qssData->lp;
  double *lqu = qssData->lqu;
  int nSD;
  const int xOrder = qssData->order;
  const int qOrder = xOrder - 1;
  const int coeffs = xOrder + 1;
  double *tq = qssTime->tq;
  double *tx = qssTime->tx;
  double *nextStateTime = qssTime->nextStateTime;
  double *x = qssData->x;
  double *q = qssData->q;
  int **SD = qssData->SD;
  const QSS_idxMap qMap = lp->qMap;
  int id = simulator->id;
  MLB_mailbox mailbox = simulator->mailbox;
  double *d = qssData->d;
  int nSZ, nLHSSt, nHD, nHZ, nLHSDsc;
  SD_eventData event = qssData->event;
  int **SZ = qssData->SZ;
  int **HD = qssData->HD;
  int **HZ = qssData->HZ;
  const QSS_idxMap eMap = lp->eMap;
  double inputTime = message.time;
  if (inputTime < simulator->previousTime)
    {
      inputTime = simulator->previousTime;
      simulator->pastEvents++;
    }
  qssTime->time = inputTime;
  qssTime->minIndex = message.index;
  qssTime->type = message.type;
  double t = qssTime->time;
  int index = qssTime->minIndex;
  int cf0 = index * coeffs, infCf0;
  QSS_StepType type = qssTime->type;
  simulator->extTrans++;
  simulator->lpTime[id] = t;
  simulator->previousTime = t;
  if (message.sendAck)
    {
      MLB_ack (mailbox, message.from, id);
    }
#ifdef DEBUG
  SD_simulationSettings settings = simulator->settings;
  SD_simulationLog simulationLog = simulator->simulationLog;
  if (settings->debug & SD_DBG_ExternalEvent)
    {
      SD_print (simulationLog, "LP %d external event: index = %d, type = %d, time = %.16lf, gvt = %.16lf, previousTime = %.16lf",id, message.index, message.type, message.time, simulator->lpTime[id], simulator->previousTime);
    }
#endif
  switch (type)
    {
    case ST_State:
      tq[index] = t;
      for (i = 0; i <= qOrder; i++)
	{
	  q[cf0 + i] = message.value[i];
	}
      // Derivative change.
      nSD = qssData->nSD[index];
      for (i = 0; i < nSD; i++)
	{
	  j = SD[index][i];
	  if (qMap[j] != NOT_ASSIGNED)
	    {
	      elapsed = t - tx[j];
	      infCf0 = j * coeffs;
	      if (elapsed > 0)
		{
		  x[infCf0] = evaluatePoly (infCf0, elapsed, x, xOrder);
		  tx[j] = t;
		}
	      FRW_recomputeDerivative (frw, qssModel, qssData, qssTime, j);
	      QA_recomputeNextTime (quantizer, j, t, nextStateTime, x, lqu, q);
	    }
	}
      nSZ = qssData->nSZ[index];
      for (i = 0; i < nSZ; i++)
	{
	  j = SZ[index][i];
	  if (eMap[j] != NOT_ASSIGNED)
	    {
	      FRW_nextEventTime (frw, qssModel, qssData, qssTime, j);
	    }
	}
      break;
    case ST_Event:
      // Handler change.
      nLHSDsc = event[index].nLHSDsc;
      for (i = 0; i < nLHSDsc; i++)
	{
	  j = event[index].LHSDsc[i];
	  d[j] = message.value[i];
	}
      nLHSSt = event[index].nLHSSt;
      for (i = 0; i < nLHSSt; i++)
	{
	  int k, h;
	  j = event[index].LHSSt[i];
	  if (qMap[j] != NOT_ASSIGNED)
	    {
	      // Trajectory change for stVars[i].
	      infCf0 = j * coeffs;
	      x[infCf0] = message.value[i * coeffs + nLHSDsc];
	      tx[j] = t;
	      qssTime->reinits->time = t;
	      qssTime->reinits->variable[qssTime->reinits->counter++] = j;
	    }
	  else if (message.value[i * coeffs + xOrder + nLHSDsc] == message.from)
	    {
	      tq[j] = t;
	      int updIdx;
	      for (updIdx = 0; updIdx <= qOrder; updIdx++)
		{
		  q[j * coeffs + updIdx] = message.value[i * coeffs + updIdx
		      + nLHSDsc];
		}
	      nSD = qssData->nSD[j];
	      for (h = 0; h < nSD; h++)
		{
		  k = qssData->SD[j][h];
		  if (qMap[k] != NOT_ASSIGNED)
		    {
		      infCf0 = k * coeffs;
		      elapsed = t - tx[k];
		      if (elapsed > 0)
			{
			  x[infCf0] = evaluatePoly (infCf0, elapsed, x, xOrder);
			  tx[k] = t;
			}
		      FRW_recomputeDerivative (frw, qssModel, qssData, qssTime,
					       k);
		      QA_recomputeNextTime (quantizer, k, t, nextStateTime, x,
					    lqu, q);
		    }
		}
	    }
	  nSZ = qssData->nSZ[j];
	  for (h = 0; h < nSZ; h++)
	    {
	      k = SZ[j][h];
	      if (eMap[k] != NOT_ASSIGNED)
		{
		  FRW_nextEventTime (frw, qssModel, qssData, qssTime, k);
		}
	    }
	}
      nHZ = qssData->nHZ[index];
      for (i = 0; i < nHZ; i++)
	{
	  j = HZ[index][i];
	  if (eMap[j] != NOT_ASSIGNED)
	    {
	      FRW_nextEventTime (frw, qssModel, qssData, qssTime, j);
	    }
	}
      nHD = qssData->nHD[index];
      for (i = 0; i < nHD; i++)
	{
	  j = HD[index][i];
	  if (qMap[j] != NOT_ASSIGNED)
	    {
	      elapsed = t - tx[j];
	      infCf0 = j * coeffs;
	      if (elapsed > 0)
		{
		  x[infCf0] = evaluatePoly (infCf0, elapsed, x, xOrder);
		  tx[j] = t;
		}
	      FRW_recomputeDerivative (frw, qssModel, qssData, qssTime, j);
	      QA_recomputeNextTime (quantizer, j, t, nextStateTime, x, lqu, q);
	    }
	}
      break;
    default:
      break;
    }
}

void
QSS_PARH_integrator (QSS_simulator simulator)
{
  int code = PAR_initLPTasks (simulator->id);
  if (code != PAR_NO_ERROR)
    {
      QSS_PAR_printParallelLog (simulator, code);
    }
  int i, j;
  double elapsed, Dt = 0, Dq = 0, Dx = 0;
  QSS_data qssData = simulator->data;
  QSS_time qssTime = simulator->time;
  FRW_framework frw = simulator->frw;
  OUT_output log = simulator->log;
  SC_scheduler scheduler = simulator->scheduler;
  QSS_model qssModel = simulator->model;
  QA_quantizer quantizer = simulator->quantizer;
  SD_output output = simulator->output;
  QSS_LP_data lp = qssData->lp;
  double t = qssTime->time;
  int index = qssTime->minIndex;
  int cf0, infCf0;
  QSS_StepType type = qssTime->type;
  int nSD, nOutputs = lp->outputs;
  QSS_dt dt = simulator->dt;
  const double ft = qssData->ft;
  const int xOrder = qssData->order;
  const int qOrder = xOrder - 1;
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
  const QSS_idxMap qMap = lp->qMap;
  int id = simulator->id;
  IBX_inbox inbox = simulator->inbox;
  MLB_mailbox mailbox = simulator->mailbox;
  double nextMessageTime;
  int synchronize = NOT_ASSIGNED;
#ifdef DEBUG
  SD_simulationSettings settings = simulator->settings;
  SD_simulationLog simulationLog = simulator->simulationLog;
#endif
  double *a = qssData->alg;
  int nSZ, nLHSSt, nRHSSt, nHD, nHZ, nLHSDsc;
  SD_eventData event = qssData->event;
  double *d = qssData->d;
  double *tmp1 = qssData->tmp1;
  int **SZ = qssData->SZ;
  int **ZS = qssData->ZS;
  int **HD = qssData->HD;
  int **HZ = qssData->HZ;
  const QSS_idxMap eMap = lp->eMap;
  t = QSS_PAR_passiveInitialization(simulator,QSS_PARH_externalEvent);
  double gvt = QSS_PAR_GVT (simulator);
  double maxAdvanceTime = gvt + QSS_dtValue (dt);
#ifdef DEBUG
  if (settings->debug & SD_DBG_StepInfo)
    {
      SD_print (simulator->simulationLog, "\nBegin Simulation:");
    }
#endif
  getTime (simulator->iTime);
  while (TRUE)
    {
      if (t == ft)
	{
	  t = QSS_PAR_passiveLP(simulator, QSS_PARH_externalEvent);
	}
      else
	{
	  // WAITFOR
	  while (t > maxAdvanceTime && gvt <= ft)
	    {
	      IBX_checkAckInbox (inbox, mailbox, id);
	      nextMessageTime = IBX_nextMessageTime (inbox);
	      if (nextMessageTime <= maxAdvanceTime)
		{
#ifdef DEBUG
		  if (settings->debug & SD_DBG_WaitFor)
		    {
		      SD_print (simulator->simulationLog, "LP %d waiting, maxAdvanceTime = % .16lf, gvt = %.16lf localTime = %.16lf",id, maxAdvanceTime, gvt, t);
		    }
#endif
		  QSS_PARH_externalEvent (simulator, IBX_nextMessage (inbox));
		  SC_update (scheduler, qssData, qssTime);
		  nextMessageTime = IBX_nextMessageTime (inbox);
		  if (nextMessageTime < qssTime->time)
		    {
		      simulator->externalEvent = TRUE;
		      qssTime->time = nextMessageTime;
		      if (nextMessageTime < simulator->previousTime)
			{
			  qssTime->time = simulator->previousTime;
			}
		    }
		  else
		    {
		      simulator->externalEvent = FALSE;
		    }
		  t = qssTime->time;
		  simulator->lpTime[id] = t;
		}
	      QSS_dtCheck (dt);
	      gvt = QSS_PAR_GVT (simulator);
	      maxAdvanceTime = gvt + QSS_dtValue (dt);
	    }
	}
      if (t >= ft)
	{
	  simulator->lpTime[id] = INF;
	  break;
	}
      if (simulator->externalEvent)
	{
	  QSS_PARH_externalEvent (simulator, IBX_nextMessage (inbox));
	  synchronize = NOT_ASSIGNED;
	}
      else
	{
#ifdef DEBUG
	  if (settings->debug & SD_DBG_StepInfo)
	    {
	      SD_print (simulationLog, "Simulation Time: %g", t);
	    }
#endif
	  index = qssTime->minIndex;
	  type = qssTime->type;
	  cf0 = index * coeffs;
	  switch (type)
	    {
	    case ST_State:
	      {
#ifdef DEBUG
		if (settings->debug & SD_DBG_StepInfo)
		  {
		    SD_print (simulationLog, "State Variable: %d", index);
		  }
#endif
		synchronize = qMap[index];
		// Internal trajectory change.
		Dt = t - tx[index];
		elapsed = x[cf0];
		advanceTime (cf0, Dt, x, xOrder);
		Dx = x[cf0] - elapsed;
		tx[index] = t;
		lqu[index] = dQRel[index] * fabs (x[cf0]);
		if (lqu[index] < dQMin[index])
		  {
		    lqu[index] = dQMin[index];
		  }
		Dq = lqu[index];
		QA_updateQuantizedState (quantizer, index, q, x, lqu);
		tq[index] = t;
		if (synchronize >= 0)
		  {
		    IBX_message msg;
		    msg.from = id;
		    msg.type = type;
		    msg.time = t;
		    msg.index = index;
		    msg.sendAck = TRUE;
		    for (i = 0; i <= qOrder; i++)
		      {
			msg.value[i] = q[cf0 + i];
		      }
		    int lpsBegin = lp->nLPS[synchronize];
		    int lpsEnd = lp->nLPS[synchronize + 1];
		    for (i = lpsBegin; i < lpsEnd; i++)
		      {
			simulator->messages++;
			MLB_send (mailbox, lp->lps[i], id, msg);
		      }
		  }
		QA_nextTime (quantizer, index, t, nextStateTime, x, lqu);
		// Derivative change.
		int inf = 0;
		nSD = qssData->nSD[index];
		for (i = 0; i < nSD; i++)
		  {
		    j = SD[index][i];
		    if (qMap[j] != NOT_ASSIGNED)
		      {
			elapsed = t - tx[j];
			infCf0 = j * coeffs;
			if (elapsed > 0)
			  {
			    x[infCf0] = evaluatePoly (infCf0, elapsed, x,
						      xOrder);
			    tx[j] = t;
			  }
			inf++;
		      }
		  }
		if (inf)
		  {
		    FRW_recomputeDerivatives (frw, qssModel, qssData, qssTime,
					      index);
		    QA_recomputeNextTimes (quantizer, nSD, qssData->SD[index],
					   t, nextStateTime, x, lqu, q);
		  }
		nSZ = qssData->nSZ[index];
		for (i = 0; i < nSZ; i++)
		  {
		    j = qssData->SZ[index][i];
		    if (eMap[j] != NOT_ASSIGNED)
		      {
			FRW_nextEventTime (frw, qssModel, qssData, qssTime, j);
		      }
		  }
		if (nOutputs)
		  {
		    if (output->nSO[index])
		      {
			OUT_write (log, qssData, qssTime, output);
		      }
		  }
	      }
	      break;
	    case ST_Event:
	      {
#ifdef DEBUG
		if (settings->debug & SD_DBG_StepInfo)
		  {
		    SD_print (simulationLog, "Event: %d", index);
		  }
#endif
		double zc[4];
		int s;
		int nZS = qssData->nZS[index];
		synchronize = eMap[index];
		IBX_message msg;
		for (i = 0; i < nZS; i++)
		  {
		    j = ZS[index][i];
		    elapsed = t - tq[j];
		    if (elapsed > 0)
		      {
			advanceTime (j * coeffs, elapsed, q, qOrder);
			tq[j] = t;
		      }
		  }
		qssModel->events->zeroCrossing (index, q, d, a, t, zc);
		// Zero crossing function change.
		s = sign (zc[0]);
		if (event[index].zcSign != s || xOrder == 1)
		  {
		    if (event[index].direction == 0
			|| event[index].direction == s)
		      {
			nRHSSt = event[index].nRHSSt;
			for (i = 0; i < nRHSSt; i++)
			  {
			    j = event[index].RHSSt[i];
			    elapsed = t - tq[j];
			    infCf0 = j * coeffs;
			    if (elapsed > 0)
			      {
				tmp1[infCf0] = evaluatePoly (infCf0, elapsed, q,
							     qOrder);
			      }
			    else
			      {
				tmp1[infCf0] = q[infCf0];
			      }
			  }
			nLHSSt = event[index].nLHSSt;
			for (i = 0; i < nLHSSt; i++)
			  {
			    j = event[index].LHSSt[i];
			    infCf0 = j * coeffs;
			    elapsed = t - tq[j];
			    if (elapsed > 0)
			      {
				tmp1[infCf0] = evaluatePoly (infCf0, elapsed, q,
							     qOrder);
			      }
			    else
			      {
				tmp1[infCf0] = q[infCf0];
			      }
			  }
			if (s >= 0)
			  {
			    qssModel->events->handlerPos (index, tmp1, d, a, t);
			  }
			else
			  {
			    qssModel->events->handlerNeg (index, tmp1, d, a, t);
			  }
			nLHSDsc = event[index].nLHSDsc;
			if (synchronize >= 0)
			  {
			    for (i = 0; i < nLHSDsc; i++)
			      {
				j = event[index].LHSDsc[i];
				msg.value[i] = d[j];
			      }
			  }
			for (i = 0; i < nLHSSt; i++)
			  {
			    j = event[index].LHSSt[i];
			    infCf0 = j * coeffs;
			    if (qMap[j] != NOT_ASSIGNED)
			      {
				x[infCf0] = tmp1[infCf0];
				tx[j] = t;
				lqu[j] = dQRel[j] * fabs (x[infCf0]);
				if (lqu[j] < dQMin[j])
				  {
				    lqu[j] = dQMin[j];
				  }
				QA_updateQuantizedState (quantizer, j, q, x,
							 lqu);
				tq[j] = t;
				if (synchronize >= 0)
				  {
				    int updIdx;
				    for (updIdx = 0; updIdx <= qOrder; updIdx++)
				      {
					msg.value[nLHSDsc + updIdx + i * coeffs] =
					    q[infCf0 + updIdx];
				      }
				    msg.value[nLHSDsc + xOrder + i * coeffs] =
					id;
				  }
			      }
			    else
			      {
				msg.value[nLHSDsc + i * coeffs] = tmp1[infCf0];
				int updIdx;
				for (updIdx = 1; updIdx <= xOrder; updIdx++)
				  {
				    msg.value[nLHSDsc + updIdx + i * coeffs] =
					NOT_ASSIGNED;
				  }
			      }
			  }
			qssModel->events->zeroCrossing (index, q, d, a, t, zc);
			event[index].zcSign = sign (zc[0]);
			for (i = 0; i < nLHSSt; i++)
			  {
			    j = event[index].LHSSt[i];
			    if (qMap[j] != NOT_ASSIGNED)
			      {
				QA_nextTime (quantizer, j, t, nextStateTime, x,
					     lqu);
				int k, h;
				int inf = 0;
				nSD = qssData->nSD[j];
				for (h = 0; h < nSD; h++)
				  {
				    k = SD[j][h];
				    if (qMap[k] != NOT_ASSIGNED)
				      {
					elapsed = t - tx[k];
					infCf0 = k * coeffs;
					if (elapsed > 0)
					  {
					    x[infCf0] = evaluatePoly (infCf0,
								      elapsed,
								      x,
								      xOrder);
					    tx[k] = t;
					  }
					inf++;
				      }
				  }
				if (inf)
				  {
				    FRW_recomputeDerivatives (frw, qssModel,
							      qssData, qssTime,
							      j);
				    QA_recomputeNextTimes (quantizer, nSD,
							   SD[j], t,
							   nextStateTime, x,
							   lqu, q);
				  }
				nSZ = qssData->nSZ[j];
				for (h = 0; h < nSZ; h++)
				  {
				    k = SZ[j][h];
				    if (k != index)
				      {
					FRW_nextEventTime (frw, qssModel,
							   qssData, qssTime, k);
				      }
				  }
			      }
			  }
			nHZ = qssData->nHZ[index];
			for (i = 0; i < nHZ; i++)
			  {
			    j = HZ[index][i];
			    if (j != index && eMap[j] != NOT_ASSIGNED)
			      {
				FRW_nextEventTime (frw, qssModel, qssData,
						   qssTime, j);
			      }
			  }
			nHD = qssData->nHD[index];
			int inf = 0;
			for (i = 0; i < nHD; i++)
			  {
			    j = HD[index][i];
			    if (qMap[j] != NOT_ASSIGNED)
			      {
				elapsed = t - tx[j];
				if (elapsed > 0)
				  {
				    infCf0 = j * coeffs;
				    x[infCf0] = evaluatePoly (infCf0, elapsed,
							      x, xOrder);
				    tx[j] = t;
				  }
				inf++;
				FRW_recomputeDerivative (frw, simulator->model,
							 qssData, qssTime, j);
			      }
			  }
			if (inf)
			  {
			    QA_recomputeNextTimes (quantizer, nHD, HD[index], t,
						   nextStateTime, x, lqu, q);
			  }
			if (nOutputs)
			  {
			    if (event[index].nLHSSt)
			      {
				OUT_write (log, qssData, qssTime, output);
			      }
			  }
			if (synchronize >= 0)
			  {
			    msg.from = id;
			    msg.type = type;
			    msg.time = t;
			    msg.index = index;
			    msg.sendAck = TRUE;
			    int lpsBegin = lp->nLPS[synchronize];
			    int lpsEnd = lp->nLPS[synchronize + 1];
			    for (i = lpsBegin; i < lpsEnd; i++)
			      {
				simulator->messages++;
				MLB_send (mailbox, lp->lps[i], id, msg);
			      }
			  }
		      }
		    else
		      {
			event[index].zcSign = sign (zc[0]);
			synchronize = NOT_ASSIGNED;
		      }
		  }
		else
		  {
		    synchronize = NOT_ASSIGNED;
		  }
		FRW_nextEventTime (frw, qssModel, qssData, qssTime, index);
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
		  }
		tx[j] = t;
		FRW_recomputeDerivative (frw, qssModel, qssData, qssTime, j);
		FRW_nextInputTime (frw, qssModel, qssData, qssTime, elapsed, j,
				   index);
		QA_recomputeNextTime (quantizer, j, t, nextStateTime, x, lqu,
				      q);
	      }
	      break;
	    default:
	      break;
	    }
	  simulator->totalSteps++;
	}
      simulator->previousTime = t;
      if (synchronize >= 0)
	{
	  if (QSS_dtLogStep (dt, Dq, Dx, Dt, synchronize))
	    {
	      gvt = QSS_PAR_GVT (simulator);
	      maxAdvanceTime = gvt + QSS_dtValue (dt);
	    }
	  QSS_PAR_synchronize (simulator, synchronize, QSS_PARH_externalEvent);
	}
      SC_update (scheduler, qssData, qssTime);
      if (qssTime->time == ft)
	{
	  IBX_checkInbox (inbox);
	}
      else
	{
	  IBX_checkAckInbox (inbox, mailbox, id);
	}
      nextMessageTime = IBX_nextMessageTime (inbox);
      if (nextMessageTime < qssTime->time)
	{
	  simulator->externalEvent = TRUE;
	  qssTime->time = nextMessageTime;
	  if (nextMessageTime < simulator->previousTime)
	    {
	      qssTime->time = simulator->previousTime;
	    }
	}
      else
	{
	  simulator->externalEvent = FALSE;
	}
      t = qssTime->time;
      simulator->lpTime[id] = t;
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
  QSS_PAR_removePendingMessages (simulator);
  QSS_PAR_controlPassiveLPS (simulator);
  QSS_PAR_waitFor (simulator);
  QSS_PAR_saveLog (simulator);
  QSS_PAR_printSimulationLog (simulator);
  PAR_cleanLPTask (id);
}

void *
QSS_PARH_runSimulation (void *sim)
{
  QSS_PAR_runSimulation(sim, QSS_PARH_integrator);
  return (NULL);
}

void
QSS_PARH_integrate (SIM_simulator simulate)
{
  QSS_simulator simulator = (QSS_simulator) simulate->state->sim;
  QSS_PAR_printParallelLog (
      simulator, PAR_createLPTasks (QSS_PARH_runSimulation, simulator));
  PAR_statistics (simulator);
}
#else
void
QSS_PARH_integrate (SIM_simulator simulate)
{
  return;
}
#endif
