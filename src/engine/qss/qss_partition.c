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

#include "qss_partition.h"

#include <stddef.h>

#include "../common/data.h"

#ifdef __linux__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <scotch/scotch.h>
#include <metis.h>

#include <common/patoh.h>
#include <common/utils.h>
#include <qss/qss_graph.h>

void
PRT_readPartition (PRT_partition partition, QSS_data data, char *name)
{
  char fileName[256];
  sprintf (fileName, "%s.part", name);
  FILE *file;
  int i;
  file = fopen (fileName, "r");
  if (file != NULL)
    {
      char * line = NULL;
      size_t len = 0;
      ssize_t read;
      int val;
      i = 0;
      while ((read = getline (&line, &len, file)) != -1)
	{
	  sscanf (line, "%d", &val);
	  partition->values[i++] = val;
	}
      fclose (file);
      if (line != NULL)
	{
	  free (line);
	}
      return;
    }
}


void
PRT_createPartitions (PRT_partition partition, QSS_data data, char *name)
{
  if (data->params->pm == SD_Manual)
    {
      PRT_readPartition(partition, data, name);
      return;
    }
  char fileName[256];
  char graphType[64] = "unweighted";
  int nparts = (data->params->lps == 0 ? 64 : data->params->lps);
  int nvtxs = data->states + data->events;
  FILE *file;
  int *xadj = NULL, *adjncy = NULL, *vwgt = NULL, *ewgt = NULL;
  int i, edges;
  SD_PartitionMethod pm = data->params->pm;
  if (GRP_readGraph (name, data, &xadj, &adjncy, &edges, 1, &vwgt, &ewgt, 0,
  NULL) == GRP_ReadError)
    {
      fprintf (stderr, "Could not read generated graph files.");
      abort ();
    }
  if (nvtxs > nparts)
    {
      switch (pm)
	{
	case SD_MetisCut:
	  {
	    idx_t ncon = 1, edgecut;
	    idx_t options[METIS_NOPTIONS];
	    METIS_SetDefaultOptions (options);
	    options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
	    METIS_PartGraphKway (&nvtxs, &ncon, xadj, adjncy, vwgt, NULL, ewgt,
				 &nparts, NULL, NULL, options, &edgecut,
				 partition->values);
	  }
	  break;
	case SD_MetisVol:
	  {
	    idx_t ncon = 1, edgecut;
	    idx_t options[METIS_NOPTIONS];
	    METIS_SetDefaultOptions (options);
	    options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
	    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;
	    METIS_PartGraphKway (&nvtxs, &ncon, xadj, adjncy, vwgt, NULL, ewgt,
				 &nparts, NULL, NULL, options, &edgecut,
				 partition->values);
	  }
	  break;
	case SD_HMetis:
	  {
	    char hgraphName[256];
	    pid_t pid;
	    if ((pid = fork ()))
	      {
		wait (NULL);
		sprintf (hgraphName, "%s.hmetis.part.%d", name, nparts);
		file = fopen (hgraphName, "r");
		if (file)
		  {
		    char * line = NULL;
		    size_t len = 0;
		    ssize_t read;
		    int index = 0;
		    while ((read = getline (&line, &len, file)) != -1)
		      {
			int part;
			sscanf (line, "%d", &part);
			partition->values[index++] = part;
		      }
		    fclose (file);
		  }
	      }
	    else
	      {
		strcpy (hgraphName, name);
		strcat (hgraphName, ".hmetis");
		file = fopen (hgraphName, "w");
		if (file)
		  {
		    int i, j;
		    if (vwgt == NULL)
		      {
			fprintf (file, "%d %d 1\n", edges, nvtxs);
		      }
		    else
		      {
			fprintf (file, "%d %d 11\n", edges, nvtxs);
		      }
		    int begin, end;
		    for (i = 0; i < edges; i++)
		      {
			begin = xadj[i];
			end = xadj[i + 1];
			fprintf (file, "%d ", ewgt[i]);
			for (j = begin; j < end; j++)
			  {
			    fprintf (file, "%d ", adjncy[j] + 1);
			  }
			fprintf (file, "\n");
		      }
		    if (vwgt != NULL)
		      {
			for (i = 0; i < nvtxs; i++)
			  {
			    fprintf (file, "%d\n", vwgt[i]);
			  }
		      }
		    fclose (file);
		  }
		char parts[10];
		sprintf (parts, "%d", nparts);
		execlp ("./khmetis", "./khmetis", hgraphName, parts, "5", "1",
			"1", "1", "0", "0", NULL);
		abort ();
	      }
	  }
	  break;
	case SD_Scotch:
	  {
	    // Run scotch partition
	    SCOTCH_Graph *graph_sc = SCOTCH_graphAlloc ();
	    SCOTCH_Strat *strat = SCOTCH_stratAlloc ();
	    SCOTCH_stratInit (strat);
	    if (SCOTCH_graphBuild (graph_sc, 0, nvtxs, xadj, xadj + 1, vwgt,
	    NULL,
				   edges, adjncy, ewgt) != 0)
	      {
		printf ("Error: Scotch Graph Build\n");
	      }
	    if (SCOTCH_graphPart (graph_sc, nparts, strat, partition->values)
		!= 0)
	      {
		printf ("Error: Scotch Graph Partition\n");
	      }
	    SCOTCH_stratExit (strat);
	    SCOTCH_graphFree (graph_sc);
	  }
	  break;
	case SD_Patoh:
	  {
	    int nconst = 1, edgecut, *partweights;
	    PaToH_Parameters args;
	    PaToH_Initialize_Parameters (&args, PATOH_CONPART,
	    PATOH_SUGPARAM_DEFAULT);
	    args._k = nparts;
	    partweights = (int *) malloc (args._k * nconst * sizeof(int));
	    PaToH_Alloc (&args, nvtxs, edges, nconst, vwgt, ewgt, xadj, adjncy);
	    PaToH_Part (&args, nvtxs, edges, nconst, 0, vwgt, ewgt, xadj,
			adjncy,
			NULL,
			partition->values, partweights, &edgecut);
	    free (partweights);
	    PaToH_Free ();
	  }
	  break;
	}
    }
  if (vwgt != NULL && ewgt != NULL)
    {
      sprintf(graphType,"weighted");
    }
  switch (pm)
    {
    case SD_MetisCut:
      sprintf (fileName, "%s-MetisCut-%s-%d.partition", name, graphType, nparts);
      break;
    case SD_MetisVol:
      sprintf (fileName, "%s-MetisVol-%s-%d.partition", name, graphType, nparts);
      break;
    case SD_HMetis:
      sprintf (fileName, "%s-HMetis-%s-%d.partition", name, graphType, nparts);
      break;
    case SD_Scotch:
      sprintf (fileName, "%s-Scoth-%s-%d.partition", name, graphType, nparts);
      break;
    case SD_Patoh:
      sprintf (fileName, "%s-Patoh-%s-%d.partition", name, graphType, nparts);
      break;
    }
  file = fopen (fileName, "w");
  if (file)
    {
      for (i = 0; i < nvtxs; i++)
	{
	  fprintf (file, "%d\n", partition->values[i]);
	}
      fclose (file);
    }
  free (xadj);
  free (adjncy);
  if (vwgt != NULL)
    {
      free (vwgt);
    }
  if (ewgt != NULL)
    {
      free (ewgt);
    }
}

PRT_partition
PRT_Partition (QSS_data data, char *name)
{
  idx_t nvtxs = data->states + data->events;
  int lps = data->params->lps, i;
  PRT_partition p = checkedMalloc (sizeof(*p));
  p->values = (idx_t*) checkedMalloc (nvtxs * sizeof(idx_t));
  cleanVector (p->values, 0, nvtxs);
  p->beginStates = 0;
  p->beginHandlers = data->states;
  p->endStates = data->states;
  p->endHandlers = data->states + data->events;
  p->nOutputs = (int*) checkedMalloc (nvtxs * sizeof(int));
  cleanVector (p->nOutputs, 0, nvtxs);
  p->outputs = (int**) checkedMalloc (nvtxs * sizeof(int*));
  for (i = 0; i < nvtxs; i++)
    {
      p->outputs[i] = (int*) checkedMalloc (lps * sizeof(int));
      cleanVector (p->outputs[i], 0, lps);
    }
  if (lps > 1)
    {
      PRT_createPartitions (p, data, name);
    }
#ifdef DEBUG
  printf("States: %d\n", data->states);
  printf("Handlers: %d\n", data->events);
  printf("Begin states: %d\n", p->beginStates);
  printf("End states: %d\n", p->endStates);
  printf("Begin handlers: %d\n", p->beginHandlers);
  printf("End handlers: %d\n\n", p->endHandlers);
#endif
  return (p);
}

void
PRT_freePartition (PRT_partition partition)
{
  int i, nvtxs = partition->endHandlers;
  free (partition->nOutputs);
  for (i = 0; i < nvtxs; i++)
    {
      free (partition->outputs[i]);
    }
  free (partition->outputs);
}

#else

PRT_partition
PRT_Partition (QSS_data data, char *name)
{
  return (NULL);
}

void
PRT_freePartition (PRT_partition partition)
{
  return;
}

#endif
