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

#include "settings.h"

#include <libconfig.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"

SD_Solver
_getSolver (const char *sol)
{
  if (!strcmp (sol, "QSS"))
    {
      return (SD_QSS);
    }
  else if (!strcmp (sol, "QSS2"))
    {
      return (SD_QSS2);
    }
  else if (!strcmp (sol, "QSS3"))
    {
      return (SD_QSS3);
    }
  else if (!strcmp (sol, "LIQSS"))
    {
      return (SD_LIQSS);
    }
  else if (!strcmp (sol, "LIQSS2"))
    {
      return (SD_LIQSS2);
    }
  else if (!strcmp (sol, "LIQSS3"))
    {
      return (SD_LIQSS3);
    }
  else if (!strcmp (sol, "CQSS"))
    {
      return (SD_CQSS);
    }
  else if (!strcmp (sol, "DASSL"))
    {
      return (SD_DASSL);
    }
  else if (!strcmp (sol, "DOPRI"))
    {
      return (SD_DOPRI);
    }
  else if (!strcmp (sol, "QSS4"))
    {
      return (SD_QSS4);
    }
  return (SD_QSS);
}

SD_PartitionMethod
_getPartitionMethod (const char *sol)
{
  if (!strcmp (sol, "MetisCut"))
    {
      return (SD_MetisCut);
    }
  else if (!strcmp (sol, "MetisVol"))
    {
      return (SD_MetisVol);
    }
  else if (!strcmp (sol, "HMetis"))
    {
      return (SD_HMetis);
    }
  else if (!strcmp (sol, "Scotch"))
    {
      return (SD_Scotch);
    }
  else if (!strcmp (sol, "Patoh"))
    {
      return (SD_Patoh);
    }
  else if (!strcmp (sol, "Manual"))
    {
      return (SD_Manual);
    }
  return (SD_MetisCut);
}

SD_DtSynch
_getDtSynch (const char *sol)
{
  if (!strcmp (sol, "SD_DT_Fixed"))
    {
      return (SD_DT_Fixed);
    }
  else if (!strcmp (sol, "SD_DT_Asynchronous"))
    {
      return (SD_DT_Adaptive);
    }
  return (SD_DT_Fixed);
}

int
_getOrder (SD_Solver sol)
{
  switch (sol)
    {
    case SD_QSS:
    case SD_CQSS:
    case SD_LIQSS:
    case SD_DASSL:
    case SD_DOPRI:
      return (1);
    case SD_QSS2:
    case SD_LIQSS2:
      return (2);
    case SD_QSS3:
    case SD_LIQSS3:
      return (3);
    case SD_QSS4:
      return (4);
    }
  return (1);
}

SET_settings
SET_Settings (char *fname)
{
  config_t cfg, *cf;
  const config_setting_t *lists;
  int count, n, ires;
  double dres;
  const char *sol;
  cf = &cfg;
  config_init (cf);
  if (!config_read_file (cf, fname))
    {
      config_destroy (cf);
      fprintf (stderr, "Initial configuration file %s not found!\n", fname);
      abort ();
      return (NULL);
    }
  SET_settings p = checkedMalloc (sizeof(*p));
  p->minstep = 0;
  p->zchyst = 0;
  p->derdelta = 0;
  p->it = 0;
  p->ft = 0;
  p->dt = 0;
  p->dqmin = NULL;
  p->dqrel = NULL;
  p->symdiff = 1;
  p->lps = 0;
  p->nodesize = 0;
  p->order = 1;
  p->solver = SD_QSS3;
  p->nDQMin = 0;
  p->nDQRel = 0;
  p->pm = SD_MetisCut;
  p->dtSynch = SD_DT_Adaptive;
  if (config_lookup_float (cf, "minstep", &dres))
    {
      if (dres == 0)
	{
	  p->minstep = MIN_STEP;
	}
      else
	{
	  p->minstep = dres;
	}
    }
  if (config_lookup_float (cf, "dt", &dres))
    {
      p->dt = dres;
    }
  if (config_lookup_float (cf, "zchyst", &dres))
    {
      if (dres == 0)
	{
	  p->zchyst = ZC_HYST;
	}
      else
	{
	  p->zchyst = dres;
	}
    }
  if (config_lookup_float (cf, "derdelta", &dres))
    {
      if (dres == 0)
	{
	  p->derdelta = DER_DELTA;
	}
      else
	{
	  p->derdelta = dres;
	}
    }
  if (config_lookup_float (cf, "it", &dres))
    {
      p->it = dres;
    }
  if (config_lookup_float (cf, "ft", &dres))
    {
      p->ft = dres;
    }
  if (config_lookup_int (cf, "symdiff", &ires))
    {
      p->symdiff = ires;
    }
  if (config_lookup_int (cf, "lps", &ires))
    {
      p->lps = ires;
    }
  if (config_lookup_int (cf, "nodesize", &ires))
    {
      if (ires == 0)
	{
	  p->nodesize = NODE_SIZE;
	}
      else
	{
	  p->nodesize = ires;
	}
    }
  if (config_lookup_string (cf, "sol", &sol))
    {
      p->solver = _getSolver (sol);
      p->order = _getOrder (p->solver);
    }
  if (config_lookup_string (cf, "partitionMethod", &sol))
    {
      p->pm = _getPartitionMethod (sol);
    }
  if (config_lookup_string (cf, "dtSynch", &sol))
    {
      p->dtSynch = _getDtSynch (sol);
    }
  lists = config_lookup (cf, "dqmin");
  count = config_setting_length (lists);
  p->dqmin = checkedMalloc (count * sizeof(double));
  for (n = 0; n < count; n++)
    {
      p->dqmin[n] = config_setting_get_float_elem (lists, n);
    }
  p->nDQMin = count;
  lists = config_lookup (cf, "dqrel");
  count = config_setting_length (lists);
  p->dqrel = checkedMalloc (count * sizeof(double));
  for (n = 0; n < count; n++)
    {
      p->dqrel[n] = config_setting_get_float_elem (lists, n);
    }
  p->nDQRel = count;
  config_destroy (cf);
  return ((p));
}

void
freeSettings (SET_settings settings)
{
  if (settings->dqmin != NULL)
    {
      free (settings->dqmin);
    }
  if (settings->dqrel != NULL)
    {
      free (settings->dqrel);
    }
  free (settings);
}
