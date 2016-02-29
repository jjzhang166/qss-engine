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


#ifdef  __linux__
#define _GNU_SOURCE
#define __USE_GNU
#include <fenv.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>

void
fpe_handler (int dummy)
{
  printf ("Floating point exception\n");
  abort ();
}
#else

#include <stdio.h>
#include <stdlib.h>
#endif

#include <common/random.h>
#include <qss/qss_commands.h>
#include <qss/qss_simulator.h>
#include <qss/qss_model.h>

QSS_simulator
QSS_Simulator ()
{
  QSS_simulator p = checkedMalloc (sizeof(*p));
  p->quantizer = NULL;
  p->log = NULL;
  p->scheduler = NULL;
  p->frw = NULL;
  p->data = NULL;
  p->time = NULL;
  p->model = NULL;
  p->output = NULL;
  p->settings = NULL;
  p->simulationLog = NULL;
  p->dt = NULL;
  p->iTime = checkedMalloc (sizeof(*(p->iTime)));
  p->sTime = checkedMalloc (sizeof(*(p->sTime)));
  p->sdTime = checkedMalloc (sizeof(*(p->sdTime)));
  p->initTime = 0;
  p->simulationTime = 0;
  p->saveTime = 0;
  p->totalSteps = 1;
  p->reinits = 0;
  p->memory = 0;
  p->sequentialMemory = 0;
  p->extTrans = 0;
  p->pastEvents = 0;
  p->id = 0;
  p->lpTime = NULL;
  p->lpDtMin = NULL;
  p->previousTime = 0;
  p->externalEvent = FALSE;
  p->mailbox = NULL;
  p->inbox = NULL;
  p->ack = NULL;
  p->stats = NULL;
  p->lps = NULL;
  p->messages = 0;
  p->messagesTime = 0;
  p->dtSteps = 0;
  p->dtSynch = NULL;
  return (p);
}

void
QSS_freeSimulator (QSS_simulator simulator)
{
  SD_freeSimulationLog (simulator->simulationLog);
  QSS_freeTime (simulator->time, simulator->data->events,
		simulator->data->inputs);
  if (simulator->settings->parallel)
    {
      if (simulator->stats != NULL)
	{
	  SD_freeOutput (simulator->output, simulator->data->states,
			 simulator->data->discretes);
	  QSS_freeModel (simulator->model);
	  QSS_freeData (simulator->data);
	  SD_freeSimulationSettings (simulator->settings);
	  free (simulator->lpTime);
	  QSS_LP_freeDataArray (simulator->lps);
	  MLB_freeMailbox (simulator->mailbox);
	  SD_freeStatistics (simulator->stats);
	}
      else
	{
	  QA_freeQuantizer (simulator->quantizer);
	  OUT_freeOutput (simulator->log);
	  SC_freeScheduler (simulator->scheduler);
	  FRW_freeFramework (simulator->frw);
	  QSS_freeDt (simulator->dt);
	}
    }
  else
    {
      QA_freeQuantizer (simulator->quantizer);
      OUT_freeOutput (simulator->log);
      SC_freeScheduler (simulator->scheduler);
      FRW_freeFramework (simulator->frw);
      SD_freeOutput (simulator->output, simulator->data->states,
		     simulator->data->discretes);
      QSS_freeData (simulator->data);
      QSS_freeModel (simulator->model);
      SD_freeSimulationSettings (simulator->settings);
    }
  free (simulator->iTime);
  free (simulator->sTime);
  free (simulator->sdTime);
  free (simulator);
}

void
QSS_simulatorEnd (SIM_simulator simulate)
{
  QSS_simulator simulator = (QSS_simulator) simulate->state->sim;
  QSS_freeSimulator (simulator);
  freeRandom ();
}

void
QSS_simulate (SIM_simulator simulate)
{
  Random ();
  QSS_simulator simulator = (QSS_simulator) simulate->state->sim;
  INT_integrator integrator = INT_Integrator (simulate);
#ifdef __linux__
  signal (SIGFPE, fpe_handler);
  feenableexcept (FE_DIVBYZERO);
#endif
  getTime (simulator->iTime);
  QSS_initializeDataStructs (simulator);
  getTime (simulator->sdTime);
  subTime (simulator->sdTime, simulator->iTime);
  simulator->initTime = getTimeValue (simulator->sdTime);
  QSS_CMD_alloc(simulator);
  INT_initialize (integrator, simulate);
  INT_integrate (integrator, simulate);
  INT_freeIntegrator (integrator);
}

void
QSS_initSimulator (SIM_simulator simulator)
{
  simulator->state->sim = (void*) QSS_Simulator ();
  ((QSS_simulator) simulator->state->sim)->settings =
      simulator->state->settings;
  simulator->ops->simulate = QSS_simulate;
  simulator->ops->freeSimulator = QSS_simulatorEnd;
}

