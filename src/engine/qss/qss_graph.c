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

#include <stdlib.h>
#include <common/utils.h>
#include <qss/qss_graph.h>

int
GRP_createGraph (QSS_data data, int **xadj, int **adjncy, int rwgt, FILE *wFile,
		 int **ewgt, int rhe, int **hevars)
{
  int states = data->states, events = data->events, nvtxs = states + events,
      edges = 0, i, j;
  if (data->params->pm == SD_HMetis || data->params->pm == SD_Patoh)
    {
      vMatrix hMatrix = VMatrix (nvtxs, 32);
      for (i = 0; i < states; i++)
	{
	  int nSD = data->nSD[i];
	  int nSZ = 0;
	  for (j = 0; j < nSD; j++)
	    {
	      int k = data->SD[i][j];
	      if (k != i)
		{
		  insert (hMatrix[i], k);
		}
	    }
	  if (data->nSZ != NULL)
	    {
	      nSZ = data->nSZ[i];
	      for (j = 0; j < nSZ; j++)
		{
		  int k = data->SZ[i][j];
		  insert (hMatrix[i], states + k);
		}
	    }
	  if (nSD > 0 || nSZ > 0)
	    {
	      edges++;
	      insert (hMatrix[i], i);
	    }
	}
      for (i = 0; i < events; i++)
	{
	  int nHD = data->nHD[i];
	  for (j = 0; j < nHD; j++)
	    {
	      int k = data->HD[i][j];
	      insert (hMatrix[states + i], k);
	    }
	  int nHZ = data->nHZ[i];
	  for (j = 0; j < nHZ; j++)
	    {
	      int k = data->HZ[i][j];
	      if (k != i)
		{
		  insert (hMatrix[states + i], states + k);
		}
	    }
	  int nLHSSt = data->event[i].nLHSSt;
	  for (j = 0; j < nLHSSt; j++)
	    {
	      int k = data->event[i].LHSSt[j];
	      insert (hMatrix[states + i], k);
	    }
	  if (nHD > 0 || nHZ > 0 || nLHSSt > 0)
	    {
	      insert (hMatrix[states + i], states + i);
	      edges++;
	    }
	  int nRHSSt = data->event[i].nRHSSt;
	  for (j = 0; j < nRHSSt; j++)
	    {
	      int k = data->event[i].RHSSt[j];
	      insert (hMatrix[k], states + i);
	      if (data->nSD[k] == 0 && data->nSZ[k] == 0)
		{
		  insert (hMatrix[k], k);
		  edges++;
		}
	    }
	}
      int pinsSize = 0;
      for (i = 0; i < nvtxs; i++)
	{
	  for (j = vectorFirst (hMatrix[i]); !vectorEnd (hMatrix[i]); j =
	      vectorNext (hMatrix[i]))
	    {
	      pinsSize++;
	    }
	}
      *xadj = (int*) checkedMalloc ((edges + 1) * sizeof(int));
      *adjncy = (int*) checkedMalloc (pinsSize * sizeof(int));
      if (rhe)
	{
	  *hevars = (int*) checkedMalloc (edges * sizeof(int));
	}
      xadj[0][0] = 0;
      int end = 0;
      if (rwgt)
	{
	  if (wFile == NULL)
	    {
	      *ewgt = (int*) checkedMalloc (edges * sizeof(int));
	      edges = 0;
	      for (i = 0; i < states; i++)
		{
		  int empty = 1;
		  for (j = vectorFirst (hMatrix[i]); !vectorEnd (hMatrix[i]);
		      j = vectorNext (hMatrix[i]))
		    {
		      adjncy[0][end++] = j;
		      empty = 0;
		    }
		  if (!empty)
		    {
		      xadj[0][edges + 1] = end;
		      ewgt[0][edges] = 1;
		      edges++;
		    }
		}
	      for (i = states; i < nvtxs; i++)
		{
		  int empty = 1;
		  for (j = vectorFirst (hMatrix[i]); !vectorEnd (hMatrix[i]);
		      j = vectorNext (hMatrix[i]))
		    {
		      adjncy[0][end++] = j;
		      empty = 0;
		    }
		  if (!empty)
		    {
		      xadj[0][edges + 1] = end;
		      ewgt[0][edges] = 10;
		      edges++;
		    }
		}
	    }
	  else
	    {
	      edges = 0;
	      for (i = 0; i < nvtxs; i++)
		{
		  int empty = 1;
		  for (j = vectorFirst (hMatrix[i]); !vectorEnd (hMatrix[i]);
		      j = vectorNext (hMatrix[i]))
		    {
		      adjncy[0][end++] = j;
		      empty = 0;
		    }
		  if (!empty)
		    {
		      xadj[0][edges + 1] = end;
		      edges++;
		    }
		}
	      *ewgt = (int *) checkedMalloc (edges * sizeof(int));
	      if (fread (ewgt[0], sizeof(int), edges, wFile) != edges)
		{
		  free (*ewgt);
		  *ewgt = NULL;
		}
	      fclose (wFile);
	    }
	}
      else
	{
	  if (rhe)
	    {
	      edges = 0;
	      for (i = 0; i < nvtxs; i++)
		{
		  int empty = 1;
		  for (j = vectorFirst (hMatrix[i]); !vectorEnd (hMatrix[i]);
		      j = vectorNext (hMatrix[i]))
		    {
		      adjncy[0][end++] = j;
		      empty = 0;
		    }
		  if (!empty)
		    {
		      xadj[0][edges + 1] = end;
		      hevars[0][edges] = i;
		      edges++;
		    }
		}
	    }
	  else
	    {
	      edges = 0;
	      for (i = 0; i < nvtxs; i++)
		{
		  int empty = 1;
		  for (j = vectorFirst (hMatrix[i]); !vectorEnd (hMatrix[i]);
		      j = vectorNext (hMatrix[i]))
		    {
		      adjncy[0][end++] = j;
		      empty = 0;
		    }
		  if (!empty)
		    {
		      xadj[0][edges + 1] = end;
		      edges++;
		    }
		}
	    }
	}
      freeVMatrix (hMatrix, nvtxs);
    }
  else
    {
      vMatrix adjMatrix = VMatrix (nvtxs, 32);
      for (i = 0; i < states; i++)
	{
	  int nSD = data->nSD[i];
	  for (j = 0; j < nSD; j++)
	    {
	      int k = data->SD[i][j];
	      if (k != i)
		{
		  insert (adjMatrix[i], k);
		  insert (adjMatrix[k], i);
		  edges += 2;
		}
	    }
	  if (data->nSZ != NULL)
	    {
	      int nSZ = data->nSZ[i];
	      for (j = 0; j < nSZ; j++)
		{
		  int k = data->SZ[i][j];
		  insert (adjMatrix[i], states + k);
		  insert (adjMatrix[states + k], i);
		  edges += 2;
		}
	    }
	}
      for (i = 0; i < events; i++)
	{
	  int nHD = data->nHD[i];
	  for (j = 0; j < nHD; j++)
	    {
	      int k = data->HD[i][j];
	      insert (adjMatrix[states + i], k);
	      insert (adjMatrix[k], states + i);
	      edges += 2;
	    }
	  int nHZ = data->nHZ[i];
	  for (j = 0; j < nHZ; j++)
	    {
	      int k = data->HZ[i][j];
	      if (k != i)
		{
		  insert (adjMatrix[states + i], states + k);
		  insert (adjMatrix[states + k], states + i);
		  edges += 2;
		}
	    }
	  int nLHSSt = data->event[i].nLHSSt;
	  for (j = 0; j < nLHSSt; j++)
	    {
	      int k = data->event[i].LHSSt[j];
	      insert (adjMatrix[states + i], k);
	      insert (adjMatrix[k], states + i);
	      edges += 2;
	    }
	  int nRHSSt = data->event[i].nRHSSt;
	  for (j = 0; j < nRHSSt; j++)
	    {
	      int k = data->event[i].RHSSt[j];
	      insert (adjMatrix[states + i], k);
	      insert (adjMatrix[k], states + i);
	      edges += 2;
	    }
	}
      *xadj = (int*) checkedMalloc ((nvtxs + 1) * sizeof(int));
      *adjncy = (int*) checkedMalloc (edges * sizeof(int));
      xadj[0][0] = 0;
      int end = 0;
      if (rwgt)
	{
	  if (wFile == NULL)
	    {
	      *ewgt = (int*) checkedMalloc (edges * sizeof(int));
	      for (i = 0; i < states; i++)
		{
		  for (j = vectorFirst (adjMatrix[i]);
		      !vectorEnd (adjMatrix[i]); j = vectorNext (adjMatrix[i]))
		    {
		      adjncy[0][end] = j;
		      ewgt[0][end++] = 1;
		    }
		  xadj[0][i + 1] = end;
		}
	      for (i = states; i < nvtxs; i++)
		{
		  for (j = vectorFirst (adjMatrix[i]);
		      !vectorEnd (adjMatrix[i]); j = vectorNext (adjMatrix[i]))
		    {
		      adjncy[0][end] = j;
		      ewgt[0][end++] = 1;
		    }
		  xadj[0][i + 1] = end;
		}
	      edges = xadj[0][nvtxs];
	    }
	  else
	    {
	      for (i = 0; i < nvtxs; i++)
		{
		  for (j = vectorFirst (adjMatrix[i]);
		      !vectorEnd (adjMatrix[i]); j = vectorNext (adjMatrix[i]))
		    {
		      adjncy[0][end++] = j;
		    }
		  xadj[0][i + 1] = end;
		}
	      edges = xadj[0][nvtxs];
	      *ewgt = (int *) checkedMalloc (edges * sizeof(int));
	      if (fread (ewgt[0], sizeof(int), edges, wFile) != edges)
		{
		  free (*ewgt);
		  *ewgt = NULL;
		}
	      fclose (wFile);
	    }
	}
      else
	{
	  for (i = 0; i < nvtxs; i++)
	    {
	      for (j = vectorFirst (adjMatrix[i]); !vectorEnd (adjMatrix[i]);
		  j = vectorNext (adjMatrix[i]))
		{
		  adjncy[0][end++] = j;
		}
	      xadj[0][i + 1] = end;
	    }
	  edges = xadj[0][nvtxs];
	}
      freeVMatrix (adjMatrix, nvtxs);
    }
  return (edges);
}

int
GRP_readGraph (char *name, QSS_data data, int **xadj, int **adjncy, int *edges,
	       int rwgt, int **vwgt, int **ewgt, int rhe, int **hevars)
{
  char fileName[256];
  FILE *file, *wFile = NULL;
  SD_PartitionMethod pm = data->params->pm;
  int nvtxs = data->states + data->events;
  if (rwgt)
    {
      sprintf (fileName, "%s.vweights", name);
      file = fopen (fileName, "rb");
      if (file)
	{
	  *vwgt = (int*) checkedMalloc (nvtxs * sizeof(int));
	  if (fread (vwgt[0], sizeof(int), nvtxs, file) != nvtxs)
	    {
	      fclose (file);
	      return (GRP_ReadError);
	    }
	  fclose (file);
	  file = NULL;
	}
    }
  if (pm == SD_Scotch || pm == SD_MetisCut || pm == SD_MetisVol)
    {
      if (rwgt)
	{
	  sprintf (fileName, "%s.eweights", name);
	  wFile = fopen (fileName, "rb");
	}
      sprintf (fileName, "%s.graph", name);
      file = fopen (fileName, "rb");
      if (file)
	{
	  *xadj = (int*) checkedMalloc ((nvtxs + 1) * sizeof(int));
	  xadj[0][0] = 0;
	  if (fread (&(xadj[0][1]), sizeof(int), nvtxs, file) != nvtxs)
	    {
	      fclose (file);
	      return (GRP_ReadError);
	    }
	  edges[0] = xadj[0][nvtxs];
	  *adjncy = (int *) checkedMalloc (edges[0] * sizeof(int));
	  if (fread (adjncy[0], sizeof(int), edges[0], file) != edges[0])
	    {
	      fclose (file);
	      return (GRP_ReadError);
	    }
	  if (wFile != NULL)
	    {
	      *ewgt = (int *) checkedMalloc (edges[0] * sizeof(int));
	      if (fread (ewgt[0], sizeof(int), edges[0], wFile) != edges[0])
		{
		  fclose (wFile);
		  return (GRP_ReadError);
		}
	      fclose (wFile);
	    }
	  fclose (file);
	}
      else
	{
	  edges[0] = GRP_createGraph (data, xadj, adjncy, rwgt, wFile, ewgt,
				      rhe, hevars);
	}
    }
  else if (pm == SD_Patoh || pm == SD_HMetis)
    {
      if (rwgt)
	{
	  sprintf (fileName, "%s.heweights", name);
	  wFile = fopen (fileName, "rb");
	}
      sprintf (fileName, "%s.hgraph", name);
      file = fopen (fileName, "rb");
      if (file)
	{
	  if (fread (edges, sizeof(int), 1, file) != 1)
	    {
	      fclose (file);
	      return (GRP_ReadError);
	    }
	  *xadj = (int*) checkedMalloc ((edges[0] + 1) * sizeof(int));
	  xadj[0][0] = 0;
	  if (fread (&(xadj[0][1]), sizeof(int), edges[0], file) != edges[0])
	    {
	      fclose (file);
	      return (GRP_ReadError);
	    }
	  int npins = xadj[0][edges[0]];
	  *adjncy = (int*) checkedMalloc (npins * sizeof(int));
	  if (fread (adjncy[0], sizeof(int), npins, file) != npins)
	    {
	      fclose (file);
	      return (GRP_ReadError);
	    }
	  if (rhe)
	    {
	      *hevars = (int*) checkedMalloc (edges[0] * sizeof(int));
	      if (fread (hevars[0], sizeof(int), edges[0], file) != edges[0])
		{
		  fclose (file);
		  return (GRP_ReadError);
		}
	    }
	  if (wFile != NULL)
	    {
	      *ewgt = (int *) checkedMalloc (edges[0] * sizeof(int));
	      if (fread (ewgt[0], sizeof(int), edges[0], wFile) != edges[0])
		{
		  fclose (wFile);
		  return (GRP_ReadError);
		}
	      fclose (wFile);
	    }
	  fclose (file);
	}
      else
	{
	  edges[0] = GRP_createGraph (data, xadj, adjncy, rwgt, wFile, ewgt,
				      rhe, hevars);
	}
    }
  return (GRP_Success);
}
