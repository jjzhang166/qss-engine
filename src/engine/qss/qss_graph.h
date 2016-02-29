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

#ifndef QSS_GRAPH_H_
#define QSS_GRAPH_H_

#include <qss/qss_data.h>

typedef enum
{
  GRP_ReadError,
  GRP_Success
} GRP_Status;

int
GRP_readGraph (char *name, QSS_data data, int **xadj, int **adjncy, int *edges, int rwgt, int **vwgt, int **ewgt, int rhe, int **hevars);

int
GRP_createGraph (QSS_data data, int **xadj, int **adjncy, int rwgt, FILE *wFile,
		 int **ewgt, int rhe, int **hevars);

#endif /* QSS_GRAPH_H_ */
