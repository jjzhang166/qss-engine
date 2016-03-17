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

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef __MACH__
#include <mach/mach_time.h>
#define CLOCK_REALTIME 0
#define CLOCK_MONOTONIC 0
#endif
#include <gsl/gsl_errno.h>
#include <gsl/gsl_poly.h>

#include <common/utils.h>

void *
checkedMalloc (unsigned long long len)
{
  void *p = malloc (len);
  if (!p)
    {
      fprintf (stderr, "\nRan out of memory!\n");
      exit (1);
    }
  return ((p));
}

#ifdef _WIN32
double
getTimeValue (struct timeval *te)
#else
double
getTimeValue (struct timespec *te)
#endif
{
#ifdef _WIN32
  return (te->tv_usec);
#else
  return ((te->tv_sec * 1e3 + te->tv_nsec / 1e6));
#endif
}

#ifdef _WIN32
void
getTime (struct timeval *te)
#else
void
getTime (struct timespec *te)
#endif
{
#ifdef __linux__
  if (clock_gettime (CLOCK_THREAD_CPUTIME_ID, te) == -1)
    {
      fprintf (stderr, "Error: Can't get CPU PROCESS time.");
    }
#elif defined (_WIN32)
  te->tv_sec = 0;
  te->tv_usec = (1000 * clock()) / CLOCKS_PER_SEC;
#elif defined (__MACH__)
  mach_timebase_info_data_t timebase;
  mach_timebase_info(&timebase);
  uint64_t time;
  time = mach_absolute_time();
  double nseconds = ((double)time * (double)timebase.numer)/((double)timebase.denom);
  double seconds = ((double)time * (double)timebase.numer)/((double)timebase.denom * 1e9);
  te->tv_sec = seconds;
  te->tv_nsec = nseconds;
#endif
  return;
}

#ifdef _WIN32
void
getTimeRes (struct timeval *te)
#else
void
getTimeRes (struct timespec *te)
#endif
{
#ifdef __linux__
  if (clock_getres (CLOCK_PROCESS_CPUTIME_ID, te) == -1)
    {
      fprintf (stderr, "Error: Can't get CPU PROCESS time.");
    }
#elif defined (_WIN32)
  te->tv_sec = 0;
  te->tv_usec = 0;
#endif
  return;
}

#ifdef _WIN32
void
subTime (struct timeval *v, struct timeval *u)
#else
void
subTime (struct timespec *v, struct timespec *u)
#endif
{
#ifdef _WIN32
  v->tv_usec -= u->tv_usec;
#else
  v->tv_sec -= u->tv_sec;
  v->tv_nsec -= u->tv_nsec;
  if (v->tv_nsec < 0)
    {
      v->tv_sec--;
      v->tv_nsec += 1000000000;
    }
#endif
}

int
sign (double x)
{
  if (x >= 0)
    {
      return (1);
    }
  else
    {
      return (-1);
    }
}

double
minPosRoot (double *coeff, int order)
{
  double mpr = -1;
  switch (order)
    {
    case 0:
      mpr = INF;
      break;
    case 1:
      if (coeff[1] == 0)
	{
	  mpr = INF;
	}
      else
	{
	  mpr = -coeff[0] / coeff[1];
	}
      ;
      if (mpr < 0)
	{
	  mpr = INF;
	}
      break;
    case 2:
      if (coeff[2] == 0 || (1000 * fabs (coeff[2])) < fabs (coeff[1]))
	{
	  if (coeff[1] == 0)
	    {
	      mpr = INF;
	    }
	  else
	    {
	      mpr = -coeff[0] / coeff[1];
	    };
	  if (mpr < 0)
	    {
	      mpr = INF;
	    }
	}
      else
	{
	  double disc;
	  disc = coeff[1] * coeff[1] - 4 * coeff[2] * coeff[0];
	  if (disc < 0)
	    {
	      //no real roots
	      mpr = INF;
	    }
	  else
	    {
	      double sd, r1;
	      sd = sqrt (disc);
	      r1 = (-coeff[1] + sd) / 2 / coeff[2];
	      if (r1 > 0)
		{
		  mpr = r1;
		}
	      else
		{
		  mpr = INF;
		};
	      r1 = (-coeff[1] - sd) / 2 / coeff[2];
	      if ((r1 > 0) && (r1 < mpr))
		{
		  mpr = r1;
		}
	    };
	}
      ;
      break;
    case 3:
      if ((coeff[3] == 0) || (1000 * fabs (coeff[3]) < fabs (coeff[2])))
	{
	  mpr = minPosRoot (coeff, 2);
	}
      else if (coeff[0] == 0)
	{
	  if (coeff[1] == 0)
	    {
	      mpr = -coeff[2] / coeff[3];
	    }
	  else
	    {
	      coeff[0] = coeff[1];
	      coeff[1] = coeff[2];
	      coeff[2] = coeff[3];
	      mpr = minPosRoot (coeff, 2);
	    }
	}
      else
	{
	  double q, r, disc;
	  q = (3 * coeff[3] * coeff[1] - coeff[2] * coeff[2]) / 9 / coeff[3]
	      / coeff[3];
	  r = (9 * coeff[3] * coeff[2] * coeff[1]
	      - 27 * coeff[3] * coeff[3] * coeff[0]
	      - 2 * coeff[2] * coeff[2] * coeff[2]) / 54 / coeff[3] / coeff[3]
	      / coeff[3];
	  disc = q * q * q + r * r;
	  mpr = INF;
	  if (disc >= 0)
	    {
	      //only one real root
	      double sd, s, t, r1;
	      sd = sqrt (disc);
	      if (r + sd > 0)
		{
		  s = pow (r + sd, 1.0 / 3);
		}
	      else
		{
		  s = -pow (fabs (r + sd), 1.0 / 3);
		};
	      if (r - sd > 0)
		{
		  t = pow (r - sd, 1.0 / 3);
		}
	      else
		{
		  t = -pow (fabs (r - sd), 1.0 / 3);
		};
	      r1 = s + t - coeff[2] / 3 / coeff[3];
	      if (r1 > 0)
		{
		  mpr = r1;
		}
	    }
	  else
	    {
	      //three real roots
	      double rho, th, rho13, costh3, sinth3, spt, smti32, r1;
	      rho = sqrt (-q * q * q);
	      th = acos (r / rho);
	      rho13 = pow (rho, 1.0 / 3);
	      costh3 = cos (th / 3);
	      sinth3 = sin (th / 3);
	      spt = rho13 * 2 * costh3;
	      smti32 = -rho13 * sinth3 * sqrt (3);
	      r1 = spt - coeff[2] / 3 / coeff[3];
	      if (r1 > 0)
		{
		  mpr = r1;
		}
	      r1 = -spt / 2 - coeff[2] / 3 / coeff[3] + smti32;
	      if ((r1 > 0) && (r1 < mpr))
		{
		  mpr = r1;
		}
	      r1 = r1 - 2 * smti32;
	      if ((r1 > 0) && (r1 < mpr))
		{
		  mpr = r1;
		}
	    };
	}
      ;
      break;
    case 4:
      {
	if ((coeff[4] == 0) || (1000 * fabs (coeff[4]) < fabs (coeff[3])))
	  {
	    mpr = minPosRoot (coeff, 3);
	  }
	else
	  {
	    double a[5];
	    double z[8];
	    int i;
	    a[0] = coeff[0];
	    a[1] = coeff[1];
	    a[2] = coeff[2];
	    a[3] = coeff[3];
	    a[4] = coeff[4];
	    gsl_set_error_handler_off ();
	    gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc (
		5);
	    mpr = INF;
	    if (gsl_poly_complex_solve (a, 5, w, z) == GSL_SUCCESS)
	      {
		for (i = 0; i < 4; i++)
		  {
		    if (z[2 * i + 1] == 0 && z[2 * i] > 0)
		      {
			if (z[2 * i] < mpr)
			  {
			    mpr = z[2 * i];
			  }
		      }
		  }
	      }
	    gsl_poly_complex_workspace_free (w);
	    gsl_set_error_handler (NULL);
	  }
      }
      break;
    }
  return (mpr);
}

void
advanceTime (int i, double dt, double *p, int order)
{
  int i0 = i;
  switch (order)
    {
    case 1:
      {
	int i1 = i0 + 1;
	p[i0] = p[i0] + dt * p[i1];
	return;
      }
    case 2:
      {
	int i1 = i0 + 1, i2 = i1 + 1;
	p[i0] = p[i0] + dt * p[i1] + dt * dt * p[i2];
	p[i1] = p[i1] + 2 * dt * p[i2];
	return;
      }
    case 3:
      {
	int i1 = i0 + 1, i2 = i1 + 1, i3 = i2 + 1;
	p[i0] = p[i0] + dt * p[i1] + dt * dt * p[i2] + dt * dt * dt * p[i3];
	p[i1] = p[i1] + 2 * dt * p[i2] + 3 * dt * dt * p[i3];
	p[i2] = p[i2] + 3 * dt * p[i3];
	return;
      }
    case 4:
      {
	int i1 = i0 + 1, i2 = i1 + 1, i3 = i2 + 1, i4 = i3 + 1;
	p[i0] = p[i0] + dt * p[i1] + dt * dt * p[i2] + dt * dt * dt * p[i3]
	    + dt * dt * dt * dt * p[i4];
	p[i1] = p[i1] + 2 * p[i2] * dt + 3 * p[i3] * dt * dt
	    + 4 * p[i4] * dt * dt * dt;
	p[i2] = p[i2] + 3 * p[i3] * dt + 6 * p[i4] * dt * dt;
	p[i3] = p[i3] + 4 * p[i4] * dt;
	return;
      }
    }
}

double
evaluatePoly (int i, double dt, double *p, int order)
{
  int i0 = i;
  switch (order)
    {
    case 0:
      return (p[i0]);
    case 1:
      {
	return (p[i0] + dt * p[i0 + 1]);
      }
    case 2:
      {
	return (p[i0] + dt * p[i0 + 1] + dt * dt * p[i0 + 2]);
      }
    case 3:
      {
	return (p[i0] + dt * p[i0 + 1] + dt * dt * p[i0 + 2]
	    + dt * dt * dt * p[i0 + 3]);
      }
    case 4:
      {
	return (p[i0] + dt * p[i0 + 1] + dt * dt * p[i0 + 2]
	    + dt * dt * dt * p[i0 + 3] + dt * dt * dt * dt * p[i0 + 4]);
      }
    }
  return (INF);
}

double
evaluateVectorPoly (double dt, double *p, int order)
{
  switch (order)
    {
    case 0:
      return (p[0]);
    case 1:
      return (p[0] + dt * p[1]);
    case 2:
      return (p[0] + dt * p[1] + dt * dt * p[2]);
    case 3:
      return (p[0] + dt * p[1] + dt * dt * p[2] + dt * dt * dt * p[3]);
    case 4:
      return (p[0] + dt * p[1] + dt * dt * p[2] + dt * dt * dt * p[3]
	  + dt * dt * dt * dt * p[4]);
    }
  return (INF);
}

node
Node (int ns, int ord)
{
  int i;
  node n = (node) malloc (sizeof(*n));
  n->val = NULL;
  if (ord > 0)
    {
      n->val = (double**) malloc (ns * sizeof(double*));
      for (i = 0; i < ns; i++)
	{
	  n->val[i] = (double*) malloc (ord * sizeof(double));
	}
    }
  n->timeVal = (double*) malloc (ns * sizeof(double));
  n->next = NULL;
  n->used = 0;
  return (n);
}

list
List (int ns, int ord)
{
  node n = Node (ns, ord);
  list l = checkedMalloc (sizeof(*l));
  l->begin = n;
  l->end = n;
  l->nodeSize = ns;
  l->order = ord;
  return (l);
}

void
append (list l, double t, double *add)
{
  int i, order = l->order, used = l->end->used;
  if (used < l->nodeSize)
    {
      l->end->timeVal[used] = t;
      for (i = 0; i < order; i++)
	{
	  l->end->val[used][i] = add[i];
	}
      l->end->used++;
    }
  else
    {
      node n = (node) malloc (sizeof(*n));
      int nodeSize = l->nodeSize;
      n->val = NULL;
      if (order > 0)
	{
	  n->val = (double**) malloc (nodeSize * sizeof(double*));
	  for (i = 0; i < nodeSize; i++)
	    {
	      n->val[i] = (double*) malloc (order * sizeof(double));
	    }
	}
      n->timeVal = (double*) malloc (nodeSize * sizeof(double));
      n->used = 1;
      n->timeVal[0] = t;
      for (i = 0; i < order; i++)
	{
	  n->val[0][i] = add[i];
	}
      n->next = NULL;
      l->end->next = n;
      l->end = n;
    }
}

/** Vector data structure */

int
vectorCompare (const void * a, const void * b)
{
  return (*(int*) a - *(int*) b);
}

vector
Vector (int cant)
{
  vector v = checkedMalloc (sizeof(*v));
  v->values = (int *) checkedMalloc (cant * sizeof(int));
  cleanVector (v->values, NOT_ASSIGNED, cant);
  v->size = cant;
  v->inc = cant;
  v->used = 0;
  v->ordered = FALSE;
  v->iter = 0;
  return (v);
}

void
freeVector (vector v)
{
  if (v->values != NULL)
    {
      free (v->values);
    }
  free (v);
}

void
insert (vector v, int value)
{
  int used = v->used;
  if (used < v->size)
    {
      v->values[used] = value;
      v->used++;
    }
  else
    {
      int size = v->size;
      v->size += v->inc;
      int *tmp = (int *) checkedMalloc (v->size * sizeof(int));
      cleanVector (tmp, NOT_ASSIGNED, v->size);
      memcpy (tmp, v->values, size * sizeof(int));
      free (v->values);
      v->values = tmp;
      v->values[used] = value;
      v->used++;
    }
  v->ordered = FALSE;
}

int
vectorFirst (vector v)
{
  if (!v->ordered)
    {
      qsort (v->values, v->used, sizeof(int), vectorCompare);
      v->ordered = TRUE;
    }
  v->iter = 0;
  return (v->values[0]);
}

int
vectorNext (vector v)
{
  int tmpIter = v->iter;
  int tmpVal = v->values[tmpIter];
  while (tmpVal == v->values[++tmpIter])
    {
      ;
    }
  v->iter = tmpIter;
  return (v->values[tmpIter]);
}

bool
vectorEnd (vector v)
{
  return (v->iter >= v->used);
}

vMatrix
VMatrix (int rows, int cols)
{
  vMatrix m = (vMatrix) checkedMalloc (rows * sizeof(*m));
  int i;
  for (i = 0; i < rows; i++)
    {
      m[i] = Vector (cols);
    }
  return (m);
}

void
freeVMatrix (vMatrix m, int rows)
{
  int i;
  for (i = 0; i < rows; i++)
    {
      freeVector (m[i]);
    }
  free (m);
}

void
cleanVector (int *vector, int value, int size)
{
  memset (vector, value, size * sizeof(int));
}

void
cleanDoubleVector (double *vector, int value, int size)
{
  memset (vector, value, size * sizeof(double));
}

#ifdef SYNC_RT
struct timeval init;
void setInitRealTime()
  {
    gettimeofday(&init,NULL);
  }

double getRealTime()
  {
    // Returns RT in microseconds
    struct timeval now,diff;
    gettimeofday(&now,NULL);
    timersub(&now, &init, &diff);
    return diff.tv_sec*1e6+diff.tv_usec;
  }

void waitUntil(double until)
  {
    double until_ms=1.0e6 * until;
    do
      {
      }while (getRealTime()<until_ms);
  }
#endif

/* Bit data structures */

BIT_vector
BIT_Vector (int bits)
{
  BIT_vector b = (BIT_vector) checkedMalloc (sizeof(*b));
  b->nwords = (bits >> 5) + 1;
  b->nbits = bits - 1;
  b->words = (word_t*) checkedMalloc (sizeof(*b->words) * b->nwords);
  b->resetMask = (word_t*) checkedMalloc (sizeof(*b->resetMask) * b->nwords);
  b->mask = 0;
  b->word = 0;
  memset (b->words, 0, sizeof(*b->words) * b->nwords);
  memset (b->resetMask, 0, sizeof(*b->resetMask) * b->nwords);
  return (b);
}

void
BIT_set (BIT_vector b, int bit)
{
  b->words[bit >> 5] |= 1 << (bit % BITS_PER_WORD);
}

void
BIT_setMask (BIT_vector b, int bit)
{
  b->resetMask[bit >> 5] |= 1 << (bit % BITS_PER_WORD);
}

void
BIT_clear (BIT_vector b, int bit)
{
  b->words[bit >> 5] &= ~(1 << (bit % BITS_PER_WORD));
}

unsigned long
BIT_isSet (BIT_vector b, int bit)
{
  return (b->words[bit >> 5] & (1 << (bit % BITS_PER_WORD)));
}

void
BIT_freeVector (BIT_vector b)
{
  if (b != NULL)
    {
      if (b->words != NULL)
	{
	  free (b->words);
	}
      if (b->resetMask != NULL)
	{
	  free (b->resetMask);
	}
      free (b);
    }
}

void
BIT_print (BIT_vector b)
{
  if (b == NULL)
    return;
  int i;
  for (i = 0; i < b->nwords; i++)
    {
      word_t w = b->words[i];
      int j;
      for (j = 0; j < BITS_PER_WORD; j++)
	{
	  printf ("%d", w & 1);
	  w >>= 1;
	}
      printf (" ");
    }
  printf ("\n");
}

BIT_matrix
BIT_Matrix (int rows, int cols)
{
  BIT_matrix b = (BIT_matrix) checkedMalloc (rows * sizeof(*b));
  int i;
  for (i = 0; i < rows; i++)
    {
      b[i] = BIT_Vector (cols);
    }
  return (b);
}

void
BIT_freeMatrix (BIT_matrix matrix, int rows)
{
  int i;
  for (i = 0; i < rows; i++)
    {
      BIT_freeVector (matrix[i]);
    }
  free (matrix);
}

int
BIT_first (BIT_vector v)
{
  int i;
  for (i = 0; i < v->nwords; i++)
    {
      word_t w = v->words[i];
      int lz = ffs (w) - 1;
      if (lz != -1)
	{
	  v->mask = 1 << lz;
	  v->word = i;
	  return (lz + BITS_PER_WORD * i);
	}
    }
  return (v->nbits + 1);
}

int
BIT_next (BIT_vector v)
{
  int i;
  for (i = v->word; i < v->nwords; i++)
    {
      word_t w = v->words[i] & ~(v->mask);
      int lz = ffs (w) - 1;
      if (lz != -1)
	{
	  v->mask |= 1 << lz;
	  v->word = i;
	  return (lz + BITS_PER_WORD * i);
	}
      v->mask = 0;
    }
  return (v->nbits + 1);
}

word_t
BIT_isAnySet (BIT_vector b)
{
  return ((b->words[0] & 0xFFFFFFFF) || (b->words[1] & 0xFFFFFFFF));
}

void
BIT_clearAll (BIT_vector b)
{
  b->words[0] &= b->resetMask[0];
  b->words[1] &= b->resetMask[1];
}

void
BIT_setAll (BIT_vector b)
{
  b->words[0] |= SET_ALL;
  b->words[1] |= SET_ALL;
}

word_t
BIT_numberOfSetBits (word_t i)
{
  i = i - ((i >> 1) & 0x55555555);
  i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
  return ((((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24);
}

word_t
BIT_setBits (BIT_vector b)
{
  return (BIT_numberOfSetBits (b->words[0]) + BIT_numberOfSetBits (b->words[1]));
}

#ifdef __linux__

/**
 * Message passing data structures.
 */

IBX_inbox
IBX_Inbox (int states, int ack)
{
  IBX_inbox p = checkedMalloc (sizeof(*p));
  if (ack == 0)
    {
      p->messages = (IBX_message*) checkedMalloc (
      MAX_LPS * sizeof(IBX_message));
      p->orderedMessages = (IBX_message*) checkedMalloc (
      BUFFER_SIZE * sizeof(IBX_message));
      int i;
      for (i = 0; i < MAX_LPS; i++)
	{
	  p->messages[i].time = INF;
	  p->messages[i].sendAck = 1;
	}
      for (i = 0; i < BUFFER_SIZE; i++)
	{
	  p->orderedMessages[i].time = INF;
	  p->orderedMessages[i].type = -1;
	  p->orderedMessages[i].sendAck = 1;
	}
    }
  else
    {
      p->messages = NULL;
      p->orderedMessages = NULL;
    }
  p->received = BIT_Vector (MAX_LPS);
  pthread_mutex_init (&(p->receivedMutex), NULL);
  p->head = 0;
  p->tail = 0;
  p->size = 0;
  if (states > 0)
    {
      p->states = (double*) checkedMalloc (states * sizeof(double));
      cleanDoubleVector (p->states, 0, states);
    }
  else
    {
      p->states = NULL;
    }
  return (p);
}

void
IBX_freeInbox (IBX_inbox inbox)
{
  if (inbox->messages != NULL)
    {
      free (inbox->messages);
    }
  if (inbox->orderedMessages != NULL)
    {
      free (inbox->orderedMessages);
    }
  if (inbox->states)
    {
      free (inbox->states);
    }
  pthread_mutex_destroy (&(inbox->receivedMutex));
  free (inbox);
}

void
IBX_add (IBX_inbox inbox, int from, IBX_message message)
{
  inbox->messages[from] = message;
  pthread_mutex_lock (&(inbox->receivedMutex));
  BIT_set (inbox->received, from);
  if (message.type == 0)
    {
      inbox->states[message.index] = message.time;
    }
  pthread_mutex_unlock (&(inbox->receivedMutex));
}

void
IBX_ack (IBX_inbox inbox, int from)
{
  pthread_mutex_lock (&(inbox->receivedMutex));
  BIT_set (inbox->received, from);
  pthread_mutex_unlock (&(inbox->receivedMutex));
}

word_t
IBX_checkMail (IBX_inbox inbox)
{
  word_t ret;
  pthread_mutex_lock (&(inbox->receivedMutex));
  ret = BIT_isAnySet (inbox->received);
  pthread_mutex_unlock (&(inbox->receivedMutex));
  return (ret);
}

word_t
IBX_ackMessages (IBX_inbox inbox)
{
  word_t ret;
  pthread_mutex_lock (&(inbox->receivedMutex));
  ret = BIT_setBits (inbox->received);
  pthread_mutex_unlock (&(inbox->receivedMutex));
  return (ret);
}

int
IBX_compare (const void * a, const void * b)
{
  if (((IBX_message*) a)->time < ((IBX_message*) b)->time)
    {
      return (-1);
    }
  if (((IBX_message*) a)->time == ((IBX_message*) b)->time)
    {
      return (0);
    }
  if (((IBX_message*) a)->time > ((IBX_message*) b)->time)
    {
      return (1);
    }
  return (0);
}

void
IBX_reset (IBX_inbox inbox)
{
  pthread_mutex_lock (&(inbox->receivedMutex));
  BIT_clearAll (inbox->received);
  pthread_mutex_unlock (&(inbox->receivedMutex));
}

void
IBX_receiveMessages (IBX_inbox inbox)
{
  pthread_mutex_lock (&(inbox->receivedMutex));
  int j;
  int tail = inbox->tail;
  int size = inbox->size;
  for (j = BIT_first (inbox->received); j < MAX_LPS;
      j = BIT_next (inbox->received))
    {
      inbox->orderedMessages[tail++] = inbox->messages[j];
      size++;
      BIT_clear (inbox->received, j);
    }
  qsort (inbox->orderedMessages, tail, sizeof(IBX_message), IBX_compare);
  inbox->size = size;
  inbox->tail = size;
  inbox->head = 0;
  pthread_mutex_unlock (&(inbox->receivedMutex));
}

void
IBX_receiveAndAckMessages (IBX_inbox inbox, MLB_mailbox mailbox, int id)
{
  pthread_mutex_lock (&(inbox->receivedMutex));
  int j;
  int tail = inbox->tail;
  int size = inbox->size;
  for (j = BIT_first (inbox->received); j < MAX_LPS;
      j = BIT_next (inbox->received))
    {
      inbox->orderedMessages[tail] = inbox->messages[j];
      if (inbox->orderedMessages[tail].type == 0)
	{
	  inbox->orderedMessages[tail].sendAck = 0;
	  MLB_ack (mailbox, inbox->orderedMessages[tail].from, id);
	}
      size++;
      tail++;
      BIT_clear (inbox->received, j);
    }
 // printf("Cola: %d\n",size);
  qsort (inbox->orderedMessages, tail, sizeof(IBX_message), IBX_compare);
  inbox->size = size;
  inbox->tail = size;
  inbox->head = 0;
  pthread_mutex_unlock (&(inbox->receivedMutex));
}

double
IBX_nextMessageTime (IBX_inbox inbox)
{
  pthread_mutex_lock (&(inbox->receivedMutex));
  int head = inbox->head;
  IBX_message msg = inbox->orderedMessages[head];
 /* while (msg.type == 0 && inbox->states[msg.index] > msg.time)
    {
      inbox->orderedMessages[head].time = INF;
      inbox->orderedMessages[head].type = -1;
      inbox->head++;
      inbox->size--;
      head++;
      msg = inbox->orderedMessages[head];
    }*/
  pthread_mutex_unlock (&(inbox->receivedMutex));
  return (inbox->orderedMessages[inbox->head].time);
}

IBX_message
IBX_nextMessage (IBX_inbox inbox)
{
  int head = inbox->head;
  IBX_message msg = inbox->orderedMessages[head];
  inbox->orderedMessages[head].time = INF;
  inbox->orderedMessages[head].type = -1;
  inbox->head++;
  inbox->size--;
  return (msg);
}

void
IBX_checkInbox (IBX_inbox inbox)
{
  word_t ret = IBX_checkMail (inbox);
  if (ret)
    {
      IBX_receiveMessages (inbox);
    }
}

void
IBX_checkAckInbox (IBX_inbox inbox, MLB_mailbox mailbox, int id)
{
  word_t ret = IBX_checkMail (inbox);
  if (ret)
    {
      IBX_receiveAndAckMessages (inbox, mailbox, id);
    }
}

void
IBX_close (IBX_inbox inbox, int dir)
{
  pthread_mutex_lock (&(inbox->receivedMutex));
  BIT_setMask (inbox->received, dir);
  pthread_mutex_unlock (&(inbox->receivedMutex));
}

MLB_mailbox
MLB_Mailbox (int lps)
{
  MLB_mailbox p = checkedMalloc (sizeof(*p));
  p->inbox = (IBX_inbox**) checkedMalloc (sizeof(IBX_inbox*) * 2);
  p->inbox[0] = (IBX_inbox*) checkedMalloc (sizeof(IBX_inbox) * lps);
  p->inbox[1] = (IBX_inbox*) checkedMalloc (sizeof(IBX_inbox) * lps);
  p->size = lps;
  int i;
  for (i = 0; i < lps; i++)
    {
      p->inbox[0][i] = NULL;
      p->inbox[1][i] = NULL;
    }
  return (p);
}

void
MLB_freeMailbox (MLB_mailbox mailbox)
{
  int i, size = mailbox->size;
  for (i = 0; i < size; i++)
    {
      if (mailbox->inbox[0][i] != NULL)
	{
	  IBX_freeInbox (mailbox->inbox[0][i]);
	}
      if (mailbox->inbox[1][i] != NULL)
	{
	  IBX_freeInbox (mailbox->inbox[1][i]);
	}
    }
  free (mailbox);
}

void
MLB_send (MLB_mailbox mailbox, int to, int from, IBX_message message)
{
  IBX_add (mailbox->inbox[MSG_EVENT][to], from, message);
}

void
MLB_ack (MLB_mailbox mailbox, int to, int from)
{
  IBX_ack (mailbox->inbox[MSG_ACK][to], from);
}

void
MLB_close (MLB_mailbox mailbox, int to, int dir)
{
  IBX_ack (mailbox->inbox[MSG_ACK][to], dir);
  IBX_close (mailbox->inbox[MSG_ACK][to], dir);
}

#else

void
MLB_freeMailbox (MLB_mailbox mailbox)
  {
    return;
  }

#endif
