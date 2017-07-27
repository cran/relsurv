#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP cmpfast(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP expc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP netfastpinter(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP netfastpinter2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP netwei(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"cmpfast",        (DL_FUNC) &cmpfast,         9},
    {"expc",           (DL_FUNC) &expc,            6},
    {"netfastpinter",  (DL_FUNC) &netfastpinter,   9},
    {"netfastpinter2", (DL_FUNC) &netfastpinter2, 10},
    {"netwei",         (DL_FUNC) &netwei,          8},
    {NULL, NULL, 0}
};

void R_init_relsurv(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
