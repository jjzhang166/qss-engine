#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <common/utils.h>
#include <qss/liqss2.h>
#define TOL 2


#ifdef QSS_PARALLEL
void
LIQSS2_PAR_init (QA_quantizer quantizer, QSS_data simData, QSS_time simTime)
#else
void
LIQSS2_init (QA_quantizer quantizer, QSS_data simData, QSS_time simTime)
#endif
{
	int states = simData->states;
  int i,j;
  quantizer->state->oldDx = (double *) malloc (states * sizeof(double));
  quantizer->state->qAux = (double *) malloc (states * sizeof(double));
  quantizer->state->ltq = (double *) malloc (states * sizeof(double));
  quantizer->state->A = (double **) malloc (states * sizeof(double*));
  quantizer->state->U0 = (double **) malloc (states * sizeof(double*));
  quantizer->state->U1 = (double **) malloc (states * sizeof(double*));
  for(i = 0; i < states; i++)
  {
   quantizer->state->A[i] = (double *)malloc(states * sizeof(double));
   quantizer->state->U0[i] = (double *)malloc(states * sizeof(double));
   quantizer->state->U1[i] = (double *)malloc(states * sizeof(double));
  }
  quantizer->state->qstate = (double *) malloc (3*states * sizeof(double));
  quantizer->state->tx = (double *) malloc (states * sizeof(double));
  quantizer->state->flag2 = (int *) malloc (states * sizeof(int));
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
      quantizer->state->ltq[i] = simData->it;
      quantizer->state->flag2[i] = 0; //this flag becomes true when a future situation ddx=0 is detected.
			quantizer->state->nTime = simTime->nextStateTime;
			quantizer->state->qstate[cf0] = simData->x[cf0];
			quantizer->state->qstate[cf0+1] = 0;//simData->x[cf0+1];
			quantizer->state->qstate[cf0+2] = 0;//simData->x[cf0+2];
			quantizer->state->tx[i] = 0;
      for (j=0; j < states; j++)
      {
   		quantizer->state->A[i][j] = 0;
   		quantizer->state->U0[i][j] = 0;
   		quantizer->state->U1[i][j] = 0;
      }
    }
#ifdef QSS_PARALLEL
  quantizer->state->qMap = simData->lp->qMap;
#endif
  quantizer->state->minStep = simData->params->minStep;
  quantizer->state->lSimTime = simTime;
  
  quantizer->state->nSD = simData->nSD;
  quantizer->state->SD = simData->SD;
  quantizer->ops->recomputeNextTimes = LIQSS2_recomputeNextTimes;
  quantizer->ops->recomputeNextTime = LIQSS2_recomputeNextTime;
  quantizer->ops->nextTime = LIQSS2_nextTime;
  quantizer->ops->updateQuantizedState = LIQSS2_updateQuantizedState;
}

#ifdef QSS_PARALLEL
void
LIQSS2_PAR_recomputeNextTime (QA_quantizer quantizer, int var, double t, double *nTime,
		       double *x, double *lqu, double *q)
#else
void
LIQSS2_recomputeNextTime (QA_quantizer quantizer, int var, double t, double *nTime,
		       double *x, double *lqu, double *q)
#endif
{
	int cf0 = var * 3, cf1 = cf0 + 1, cf2 = cf1 + 1;
	double diffxq[3];
	double dt1, nT_prev;
	double **A = quantizer->state->A;	
	double **U0 = quantizer->state->U0;	
	double **U1 = quantizer->state->U1;
	int stind  = quantizer->state->lSimTime->minIndex;
	double diffQ;
	int i=var;
	bool *change = quantizer->state->change;
	int *flag2 = quantizer->state->flag2;

	if(t>0)
	{
		diffQ = q[3*stind] - quantizer->state->qAux[stind];
		if (diffQ)
		{	
			A[var][stind] = (x[cf1] - quantizer->state->oldDx[i]) / diffQ;		
		}
		U0[var][var] = x[cf1] - q[cf0]  * A[var][var];
  	U1[var][var] = 2 * x[cf2] - q[cf1] * A[var][var];
	}
	else
	{
		U0[var][var] = x[cf1];
  	U1[var][var] = 2 * x[cf2];
	}
	
	nT_prev = nTime[var];
	if(t>0)
	{
		nTime[var] = INF;
		if(flag2[var] == 1) flag2[var] = 0;
		else
		{
			if((q[cf1]-x[cf1])*x[cf2]>0)
			{
				nTime[var] = t + (q[cf1]-x[cf1])/(2*x[cf2]);
			}
		}
		diffxq[1] = q[cf1]-x[cf1];
		diffxq[2] = -x[cf2];
					 	
		if(nTime[var] < nT_prev) nTime[var] = nT_prev;	//<
 	
		if(q[cf0]!=x[cf0])
		{
			if(q[cf0] > x[cf0])
			{
				diffxq[0] = q[cf0] - x[cf0] + lqu[var]/10;
			}
			else
			{
				diffxq[0] = q[cf0] - x[cf0] - lqu[var]/10;
			}
			dt1 = t + minPosRoot(diffxq,2);
			if (dt1 < nTime[var])
			{ 
				nTime[var]=dt1;
			}
			if(q[cf0] > x[cf0])
			{
				diffxq[0] = q[cf0] - x[cf0] - TOL*lqu[var];
			}
			else
			{
				diffxq[0] = q[cf0] - x[cf0] + TOL*lqu[var];
			}
			dt1 = t + minPosRoot(diffxq,2);
			if (dt1 < nTime[var])
			{ 
				nTime[var]=dt1;
			}
		}
		else
		{
			diffxq[0] = q[cf0] - x[cf0] - 1*lqu[var];
			dt1 = t + minPosRoot(diffxq,2);
			if (dt1 < nTime[var])
			{ 
				nTime[var] = dt1;
			}
			diffxq[0] = q[cf0] - x[cf0] + 1*lqu[var];
			dt1 = t + minPosRoot(diffxq,2);
			if (dt1 < nTime[var])
			{
				nTime[var] = dt1;
			}
		}
	}
	double err1 = q[cf0] - x[cf0] + diffxq[1] * (nTime[var] - t) / 2 + diffxq[2] * pow ((nTime[var] - t) / 2, 2);
  if (fabs (err1) > 4 * fabs (lqu[var])) //3
  {
		nTime[var] = t + quantizer->state->finTime * quantizer->state->minStep;
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
      if (map[inf[i]] != NOT_ASSIGNED)
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
	int cf0 = var * 3, cf1 = cf0 + 1, cf2 = cf1 + 1;
	double diffxq[3];
	double dt1;
	double *q = quantizer->state->qstate;
	nTime[var] = INF;
	
	diffxq[1] = q[cf1] - x[cf1];
	diffxq[2] = -x[cf2];
	if(q[cf0] != x[cf0])
	{
		diffxq[0] = q[cf0] - x[cf0] - TOL*lqu[var];
 		dt1 = t + minPosRoot(diffxq,2);
 		if (dt1 < nTime[var])
 		{ 
			nTime[var] = dt1;
		}
		diffxq[0] = q[cf0] - x[cf0] + TOL*lqu[var];
	  dt1 = t + minPosRoot(diffxq,2);
	  if (dt1 < nTime[var])
	  {
			nTime[var] = dt1;
		}
	}
	else
	{
		diffxq[0] = q[cf0] - x[cf0] - 1*lqu[var];
 		dt1 = t + minPosRoot(diffxq,2);
 		if (dt1 < nTime[var])
 		{
			nTime[var] = dt1;
		}
		diffxq[0] = q[cf0] - x[cf0] + 1*lqu[var];
 		dt1 = t + minPosRoot(diffxq,2);
		if (dt1 < nTime[var])
 		{
			nTime[var] = dt1;
		}
	}
}


#ifdef QSS_PARALLEL
void
LIQSS2_PAR_solve_single(QA_quantizer quantizer, int i, double *x, double *q, double *lqu)
#else
void
LIQSS2_solve_single(QA_quantizer quantizer, int i, double *x, double *q, double *lqu)
#endif
{
	double **U0 = quantizer->state->U0;
	double **U1 = quantizer->state->U1;
	double **A = quantizer->state->A;
	int i0 = 3*i, i1 = i0 + 1, i2 = i1 + 1;
	double h;
	h=INF;
	if(A[i][i]!=0)
	{
		if(A[i][i]*(A[i][i]*x[i0]+U0[i][i])+U1[i][i]<0)
		{
			q[i0]=x[i0]+lqu[i];
		}
		else
		{
			q[i0]=x[i0]-lqu[i];
		}
		double a,b,c;
		double resolv;
		double num,den;
		a=0.5;
		b=-A[i][i]*x[i0]-U0[i][i];
		c=(q[i0]-x[i0])*U1[i][i]+0.5*(A[i][i]*q[i0]+U0[i][i])*(A[i][i]*q[i0]+U0[i][i]);
		resolv=b*b-4*a*c;
		if(resolv>0)
		{
			q[i1]=(-b-sqrt(resolv))/2/a;
			num=q[i1]-A[i][i]*q[i0]-U0[i][i];
			den=A[i][i]*q[i1]+U1[i][i];
			if(num*den>0)
			{
				h=1;
			}
			else
			{
				q[i1]=(-b+sqrt(resolv))/2/a;
				num=q[i1]-A[i][i]*q[i0]-U0[i][i];
				den=A[i][i]*q[i1]+U1[i][i];
				if(num*den>0)
				{
					h=1;
				}
			}
		}    
		if(h==INF)
		{
			if(A[i][i]*(A[i][i]*x[i0]+U0[i][i])+U1[i][i]<0)
			{
				q[i0]=x[i0]-lqu[i];
			}
			else
			{
				q[i0]=x[i0]+lqu[i];
			}
			a=0.5;
			b=-A[i][i]*x[i0]-U0[i][i];
			c=(q[i0]-x[i0])*U1[i][i]+0.5*(A[i][i]*q[i0]+U0[i][i])*(A[i][i]*q[i0]+U0[i][i]);
			resolv=b*b-4*a*c;
			if(resolv>0)
			{
				q[i1]=(-b-sqrt(resolv))/2/a;
				num=q[i1]-A[i][i]*q[i0]-U0[i][i];
				den=A[i][i]*q[i1]+U1[i][i];
				if(num*den>0)
				{
					h=1;
				}
				else
				{
					q[i1]=(-b+sqrt(resolv))/2/a;
					num=q[i1]-A[i][i]*q[i0]-U0[i][i];
					den=A[i][i]*q[i1]+U1[i][i];
					if(num*den>0)
					{
						h=1;
					}
				}
			}
		}
		if(h==INF)
		{
			q[i1]=-U1[i][i]/A[i][i];
			q[i0]=(q[i1]-U0[i][i])/A[i][i];
		}
	}
	else
	{ // A[i][i] =0
		if(x[i2]!=0)
		{
			if(x[i2]>0)
			{
				q[i0]=x[i0]-lqu[i];
			}
			else
			{
				q[i0]=x[i0]+lqu[i];
			}
			h=sqrt(2*(x[i0]-q[i0])/x[i2]);
			q[i1]=x[i1]+x[i2]*h;
		}
		else
		{
				q[i0]=x[i0];
				q[i1]=x[i1];
		}
	} 
}

#ifdef QSS_PARALLEL
void
PAR_LIQSS2_old_dx(QA_quantizer quantizer, int i, double t, int nSD, double  *x, double *tx)
#else
void
LIQSS2_old_dx(QA_quantizer quantizer, int i, double t, int nSD, double  *x, double *tx)
#endif
{
	int m, j, j0,j1,j2;
	if(t>0)
	{
		for (m = 0; m < nSD; m++)
		{
			j = quantizer->state->SD[i][m];
			j0 = j * 3;
			j1 = j0 + 1;
			j2 = j1 + 1;
			// guardado derivada anterior de xj
			quantizer->state->oldDx[j] = x[j1] + (t - tx[j]) * x[j2] * 2;
			tx[j] = t;
		}
	}
}


#ifdef QSS_PARALLEL
void
LIQSS2_PAR_updateQuantizedState (QA_quantizer quantizer, int i, double *q, double *x,
			  double *lqu)
#else
void
LIQSS2_updateQuantizedState (QA_quantizer quantizer, int i, double *q, double *x,
			  double *lqu)
#endif
{
	double t = quantizer->state->lSimTime->time;
	double **U0 = quantizer->state->U0;
	double **U1 = quantizer->state->U1;
	double **A = quantizer->state->A;
	double *tx = quantizer->state->tx;
	int i0 = i * 3, i1 = i0 + 1, i2 = i1 + 1;		
	double qi_acum;
	double elapsed;
	int m,j,j0,j1,j2;
	int nSD = quantizer->state->nSD[i];
		
  elapsed = t - quantizer->state->lSimTime->tq[i];
  quantizer->state->qAux[i] = q[i0] + elapsed * q[i1];
	// guardado derivada anterior (se hace en el for)
	quantizer->state->oldDx[i] = x[i1];
	// adelantado del resto U0ii
  elapsed = t - tx[i];
  U0[i][i] = U0[i][i] + elapsed * U1[i][i];
	tx[i] = t;
 	qi_acum = quantizer->state->qAux[i];		
	double qj_acum, ddxi;	
	double qiold[2],qjold[2];
	// guarado qi anterior	
	qiold[0] = q[i0];
	qiold[1] = q[i1];
#ifdef QSS_PARALLEL
        // solver una variable
	LIQSS2_PAR_solve_single(quantizer, i, x, q, lqu);
	// guardado derivadas viejas de todas las variables influenciadas
	LIQSS2_PAR_old_dx(quantizer,i, t, nSD, x, tx);
#else
        // solver una variable
	LIQSS2_solve_single(quantizer, i, x, q, lqu);
	// guardado derivadas viejas de todas las variables influenciadas
	LIQSS2_old_dx(quantizer,i, t, nSD, x, tx);
#endif
	// guarado de estado q para usar en nextime
	quantizer->state->qstate = q;
}
