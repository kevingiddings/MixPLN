#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _DangLoop_EM_loop(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _DangLoop_rcpparma_bothproducts(SEXP);
extern SEXP _DangLoop_rcpparma_hello_world();
extern SEXP _DangLoop_rcpparma_innerproduct(SEXP);
extern SEXP _DangLoop_rcpparma_outerproduct(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_DangLoop_EM_loop",               (DL_FUNC) &_DangLoop_EM_loop,               16},
    {"_DangLoop_rcpparma_bothproducts", (DL_FUNC) &_DangLoop_rcpparma_bothproducts,  1},
    {"_DangLoop_rcpparma_hello_world",  (DL_FUNC) &_DangLoop_rcpparma_hello_world,   0},
    {"_DangLoop_rcpparma_innerproduct", (DL_FUNC) &_DangLoop_rcpparma_innerproduct,  1},
    {"_DangLoop_rcpparma_outerproduct", (DL_FUNC) &_DangLoop_rcpparma_outerproduct,  1},
    {NULL, NULL, 0}
};

void R_init_DangLoop(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
