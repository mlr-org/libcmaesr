#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern SEXP c_cmaes_wrap_single(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP c_cmaes_wrap_batch(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"c_cmaes_wrap_single", (DL_FUNC) &c_cmaes_wrap_single, 5},
    {"c_cmaes_wrap_batch", (DL_FUNC) &c_cmaes_wrap_batch, 5},
    {NULL, NULL, 0}
};

void R_init_libcmaesr(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
