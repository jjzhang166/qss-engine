#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{ 
	if (argc != 3)
		{
			printf("./read-binary-graph NAME N\n");
			printf("where:\n");
			printf(" -- NAME is the name of the with the graph and weights files in binary format (without extension)\n");
			printf(" -- N is the number of vertices\n");
			return (0);
		}
  FILE *file, *eFile, *heFile, *vFile, *teFile, *theFile, *tvFile; 
	char *name = argv[1];
	int nvtxs;
	sscanf (argv[2],"%d",&nvtxs);
	char fileName[256];
	sprintf(fileName, "%s.graph",name);
  file = fopen (fileName, "rb");
	sprintf(fileName, "%s.vweights",name);
  vFile = fopen (fileName, "rb");
	sprintf(fileName, "%s.ewgts",name);
  eFile = fopen (fileName, "rb");
	sprintf(fileName, "%s.hewgts",name);
  heFile = fopen (fileName, "rb");
	sprintf(fileName, "%s-text.vweights",name);
  tvFile = fopen (fileName, "w");
	sprintf(fileName, "%s-text.ewgts",name);
  teFile = fopen (fileName, "w");
	sprintf(fileName, "%s-text.hewgts",name);
  theFile = fopen (fileName, "w");
  int *vwgt = (int*) malloc (nvtxs * sizeof(int));
  fread (vwgt, sizeof(int), nvtxs, vFile);
  int *xadj = (int*) malloc ((nvtxs + 1) * sizeof(int));
  xadj[0] = 0;
  fread (&(xadj[1]), sizeof(int), nvtxs, file); 
  int edges = xadj[nvtxs];
  int *ewgt = (int *) malloc (edges * sizeof(int));
  fread (ewgt, sizeof(int), edges, eFile) ;
  int i;
  for (i = 0; i < nvtxs; i++)
  	{
    	fprintf(tvFile,"%d\n",vwgt[i]);
   	}
 	for (i = 0; i < edges; i++)
  	{
    	fprintf(teFile,"%d\n",ewgt[i]);
   	}    
  free (ewgt);
  fclose (file);
	sprintf (fileName, ".hgraph");
	file = fopen (fileName, "rb");
  if (file)
 	  {
		  fread (&edges, sizeof(int), 1, file);
			ewgt = (int *) malloc (edges * sizeof(int));
		  fread (ewgt, sizeof(int), edges, heFile);
 			for (i = 0; i < edges; i++)
  			{
    			fprintf(theFile,"%d\n",ewgt[i]);
   			}
		}
 	free (vwgt);
  free (xadj);
  free (ewgt);
  fclose (file);
  fclose (eFile);
  fclose (heFile);
  fclose (vFile);
  fclose (teFile);
  fclose (theFile);
  fclose (tvFile);
  return (0);
}
