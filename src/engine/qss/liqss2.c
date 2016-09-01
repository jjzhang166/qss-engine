#include <qss/liqss2.h>

#include <math.h>
#include <stdlib.h>

#include "../common/data.h"
#include "../common/utils.h"
#include "qss_data.h"
#include "qss_quantizer.h"

#ifdef QSS_PARALLEL
void
LIQSS2_PAR_init (QA_quantizer quantizer, QSS_data simData, QSS_time simTime)
#else
void
LIQSS2_init (QA_quantizer quantizer, QSS_data simData, QSS_time simTime)
#endif
{
  int states = simData->states;
  int i;
  quantizer->state->dq = (double *) malloc (states * sizeof(double));
  quantizer->state->oldDx = (double *) malloc (states * sizeof(double));
  quantizer->state->qAux = (double *) malloc (states * sizeof(double));
  quantizer->state->a = (double *) malloc (states * sizeof(double));
  quantizer->state->u0 = (double *) malloc (states * sizeof(double));
  quantizer->state->u1 = (double *) malloc (states * sizeof(double));
  quantizer->state->lt = (double *) malloc (states * sizeof(double));
  quantizer->state->ltq = (double *) malloc (states * sizeof(double));
  quantizer->state->lquOld = (double *) malloc (states * sizeof(double));
  quantizer->state->flag2 = (int *) malloc (states * sizeof(int));
  quantizer->state->flag3 = (int *) malloc (states * sizeof(int));
  quantizer->state->flag4 = (int *) malloc (states * sizeof(int));
  quantizer->state->finTime = simData->ft;
  for (i = 0; i < states; i++)
    {
      int cf0 = i * 3;
      simData->x[cf0 + 2] = 0;
      simData->q[cf0] = simData->x[cf0];
      simData->q[cf0 + 1] = 0;
      simData->tmp1[cf0] = simData->x[cf0];
      quantizer->state->qAux[i] = simData->x[cf0];
      quantizer->state->oldDx[i] = 0;
      quantizer->state->a[i] = 0;
      quantizer->state->u0[i] = 0;
      quantizer->state->u1[i] = 0;
      quantizer->state->dq[i] = 0;
      quantizer->state->lt[i] = simData->it;
      quantizer->state->ltq[i] = simData->it;
      quantizer->state->lquOld[i] = simData->lqu[i];
      quantizer->state->flag2[i] = 0; //this flag becomes true when a future situation ddx=0 is detected.
      quantizer->state->flag3[i] = 0; //this flag becomes true after trying to provoke ddx=0.
      quantizer->state->flag4[i] = 0; //this flag becomes true after detecting a sign change in ddx.
    }
  quantizer->state->minStep = simData->params->minStep;
  quantizer->state->lSimTime = simTime;
#ifdef QSS_PARALLEL
  quantizer->state->qMap = simData->lp->qMap;
  quantizer->ops->recomputeNextTimes = LIQSS2_PAR_recomputeNextTimes;
  quantizer->ops->recomputeNextTime = LIQSS2_PAR_recomputeNextTime;
  quantizer->ops->nextTime = LIQSS2_PAR_nextTime;
  quantizer->ops->updateQuantizedState = LIQSS2_PAR_updateQuantizedState;
#else
  quantizer->ops->recomputeNextTimes = LIQSS2_recomputeNextTimes;
  quantizer->ops->recomputeNextTime = LIQSS2_recomputeNextTime;
  quantizer->ops->nextTime = LIQSS2_nextTime;
  quantizer->ops->updateQuantizedState = LIQSS2_updateQuantizedState;
#endif
}

#ifdef QSS_PARALLEL
void
LIQSS2_PAR_recomputeNextTime (QA_quantizer quantizer, int var, double t,
			  double *nTime, double *x, double *lqu, double *q)
#else
void
LIQSS2_recomputeNextTime (QA_quantizer quantizer, int var, double t,
			  double *nTime, double *x, double *lqu, double *q)
#endif
{
  int cf0 = var * 3, cf1 = cf0 + 1, cf2 = cf1 + 1;
  double *u0 = quantizer->state->u0;
  double *u1 = quantizer->state->u1;
  double *a = quantizer->state->a;
  double *dq = quantizer->state->dq;
  if (quantizer->state->ltq[var] == t)
    {
      double diffQ;
      diffQ = q[cf0] - quantizer->state->qAux[var];
      if (fabs (diffQ) > lqu[var] * 1e-6)
	{
	  a[var] = (x[cf1] - quantizer->state->oldDx[var]) / diffQ;
	  if (a[var] > 0)
	    {
	      a[var] = 0;
	    }
	}
    }
  else
    {
      quantizer->state->flag3[var] = 0;
    }
  u0[var] = x[cf1] - q[cf0] * a[var];
  u1[var] = 2 * x[cf2] - q[cf1] * a[var];
  quantizer->state->lt[var] = t;
  if (quantizer->state->flag4[var])
    {
      nTime[var] = t;
    }
  else
    {
      double diffxq[3];
      double timeaux;
      diffxq[1] = q[cf1] - x[cf1];
      diffxq[2] = -x[cf2];
      diffxq[0] = q[cf0] - dq[var] + lqu[var] - x[cf0];
      nTime[var] = t + minPosRoot (diffxq, 2);
      diffxq[0] = q[cf0] - dq[var] - lqu[var] - x[cf0];
      timeaux = t + minPosRoot (diffxq, 2);
      if (timeaux < nTime[var])
	{
	  nTime[var] = timeaux;
	}
      if (a[var] != 0 && (fabs (x[cf2]) > 1e-10)
	  && !quantizer->state->flag3[var] && !quantizer->state->flag2[var])
	{
	  double coeff[2];
	  coeff[0] = a[var] * a[var] * q[cf0] + a[var] * u0[var]
	      + u1[var];
	  coeff[1] = a[var] * a[var] * q[cf1] + a[var] * u1[var];
	  timeaux = t + minPosRoot (coeff, 1);
	  if (timeaux < nTime[var])
	    {
	      quantizer->state->flag2[var] = 1;
	      nTime[var] = timeaux;
	      quantizer->state->lquOld[var] = lqu[var];
	    }
	}
      else
	{
	  quantizer->state->flag2[var] = 0;
	}
      if (nTime[var] > quantizer->state->finTime)
	{
	  nTime[var] = quantizer->state->finTime;
	}
      double err1 = q[cf0] - x[cf0] + diffxq[1] * (nTime[var] - t) / 2
	  + diffxq[2] * pow ((nTime[var] - t) / 2, 2);
      if (fabs (err1) > 3 * fabs (lqu[var]))
	{
	  nTime[var] = t
	      + quantizer->state->finTime * quantizer->state->minStep;
	}
    }
 //           printf("time=%g: q[0]=%g, q[1]=%g, next time=%g\n",t,q[cf0],q[cf1],nTime[var]);
    if (q[cf0]*q[cf1]<0&&fabs(q[cf0])>10*lqu[var]) 
    {
        double timeaux2=-q[cf0]/q[cf1]-2*fabs(lqu[var]/q[cf1]);
        if (nTime[var]>t+timeaux2)
        {
            nTime[var]=t+timeaux2;
//           printf("time=%g: q[0]=%g, q[1]=%g, dQ=%g, new next time=%g\n",t,q[cf0],q[cf1],lqu[var],nTime[var]);
        }
    }
}

#ifdef QSS_PARALLEL
void
LIQSS2_PAR_recomputeNextTimes (QA_quantizer quantizer, int vars, int *inf, double t,
			   double *nTime, double *x, double *lqu, double *q)
#else
void
LIQSS2_recomputeNextTimes (QA_quantizer quantizer, int vars, int *inf, double t,
			   double *nTime, double *x, double *lqu, double *q)
#endif
{
  int i;
#ifdef QSS_PARALLEL
  int *map = quantizer->state->qMap;
#endif
  for (i = 0; i < vars; i++)
    {
#ifdef QSS_PARALLEL
      if (map[inf[i]] > NOT_ASSIGNED)
	{
#endif
      LIQSS2_recomputeNextTime (quantizer, inf[i], t, nTime, x, lqu, q);
#ifdef QSS_PARALLEL
    }
#endif
    }
}

#ifdef QSS_PARALLEL
void
LIQSS2_PAR_nextTime (QA_quantizer quantizer, int var, double t, double *nTime,
		 double *x, double *lqu)
#else
void
LIQSS2_nextTime (QA_quantizer quantizer, int var, double t, double *nTime,
		 double *x, double *lqu)
#endif
{
  int cf2 = var * 3 + 2;
  if (x[cf2] == 0)
    {
      nTime[var] = INF;
    }
  else
    {
      nTime[var] = t + sqrt (fabs (lqu[var] / x[cf2]));
    }
}

#ifdef QSS_PARALLEL
void
LIQSS2_PAR_updateQuantizedState (QA_quantizer quantizer, int var, double *q,
			     double *x, double *lqu)
#else
void
LIQSS2_updateQuantizedState (QA_quantizer quantizer, int var, double *q,
			     double *x, double *lqu)
#endif
{
  int cf0 = var * 3, cf1 = cf0 + 1, cf2 = cf1 + 1;
  double dx, elapsed;
  double *u0 = quantizer->state->u0;
  double *u1 = quantizer->state->u1;
  double *a = quantizer->state->a;
  double *dq = quantizer->state->dq;
  quantizer->state->flag3[var] = 0;
  elapsed = quantizer->state->lSimTime->time
      - quantizer->state->lSimTime->tq[var];
  quantizer->state->qAux[var] = q[cf0] + elapsed * q[cf1];
  quantizer->state->oldDx[var] = x[cf1];
  elapsed = quantizer->state->lSimTime->time - quantizer->state->lt[var];
  quantizer->state->ltq[var] = quantizer->state->lSimTime->time;
  u0[var] = u0[var] + elapsed * u1[var];
  if (quantizer->state->flag2[var])
    {
      lqu[var] = quantizer->state->lquOld[var];
      quantizer->state->flag2[var] = 0;
      q[cf0] = quantizer->state->qAux[var];
      if  (fabs(q[cf0] - x[cf0])>2*lqu[var]) q[cf0]=x[cf0];
    }
  else
    {
      q[cf0] = x[cf0];
    }
  if (a[var] < -1e-30)
    {
      if (x[cf2] < 0)
	{
	  dx = a[var] * a[var] * (q[cf0] + lqu[var]) + a[var] * u0[var]
	      + u1[var];
	  if (dx <= 0)
	    {
	      dq[var] = lqu[var];
	    }
	  else
	    {
	      dq[var] = (-u1[var] / a[var] / a[var]) - (u0[var] / a[var])
		  - q[cf0];
	      quantizer->state->flag3[var] = 1;
	      if (fabs (dq[var]) > lqu[var])
		{
		  dq[var] = lqu[var];
		}
	    }
	}
      else
	{
	  dx = a[var] * a[var] * (q[cf0] - lqu[var]) + a[var] * u0[var]
	      + u1[var];
	  if (dx >= 0)
	    {
	      dq[var] = -lqu[var];
	    }
	  else
	    {
	      dq[var] = -u1[var] / a[var] / a[var] - u0[var] / a[var] - q[cf0];
	      quantizer->state->flag3[var] = 1;
	      if (fabs (dq[var]) > lqu[var])
		{
		  dq[var] = -lqu[var];
		}
	    }
	}
      if (q[cf1] * x[cf1] < 0 && !quantizer->state->flag4[var]
	  && !quantizer->state->flag2[var] && !quantizer->state->flag3[var])
	{
	  if (q[cf1] < 0)
	    {
	      dq[var] = quantizer->state->qAux[var] - q[cf0]
		  - fabs (quantizer->state->lquOld[var]) * 0.1;
	    }
	  else
	    {
	      dq[var] = quantizer->state->qAux[var] - q[cf0]
		  + fabs (quantizer->state->lquOld[var]) * 0.1;
	    }
	  quantizer->state->flag4[var] = 1;
	}
      else if (quantizer->state->flag4[var])
	{
	  quantizer->state->flag4[var] = 0;
	  if (fabs (-u1[var] / a[var] / a[var] - u0[var] / a[var] - q[cf0])
	      < 3 * lqu[var])
	    {
	      dq[var] = -u1[var] / a[var] / a[var] - u0[var] / a[var] - q[cf0];
	      quantizer->state->flag3[var] = 1;
	    }
	}
    }
  else
    {
      quantizer->state->flag4[var] = 0;
      if (x[cf2] < 0)
	{
	  dq[var] = -lqu[var];
	}
      else
	{
	  dq[var] = lqu[var];
	}
    }
  if (fabs (dq[var]) > 2* lqu[var])
    {
      if (dq[var] > 0)
	{
	  dq[var] = lqu[var];
	}
      else
	{
            dq[var] = -lqu[var];
	}
    }
  q[cf0] = q[cf0] + dq[var];
  if (quantizer->state->flag3[var])
    {
      q[cf1] = a[var] * q[cf0] + u0[var];
    }
  else
    {
      q[cf1] = x[cf1];
    }
}
