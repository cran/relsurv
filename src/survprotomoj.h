/*
** $Id $
**
** Prototypes of all the survival functions
**  Including this in each routine helps prevent mismatched argument errors
*/
#include "R.h"
#include "Rinternals.h"
#define ALLOC(a,b)  R_alloc(a,b)

double **dmatrix(double *array, int ncol, int nrow);


SEXP expc(SEXP   efac2,   SEXP edims2,
		  	      SEXP   ecut2,     SEXP   expect2,
		  	      SEXP   x2, 	SEXP   y2);


SEXP netwei(   SEXP   efac2,   SEXP edims2,
	      SEXP   ecut2,     SEXP   expect2,
	      SEXP   x2, 	SEXP   y2,  SEXP status2,    SEXP times2) ;


double pystep(int nc,        int  *index,  int  *index2,   double *wt,
	      double *data,  Sint *fac,    Sint *dims,     double **cuts,
	      double step,   int  edge);

