#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL

extern SEXP c_cmaes_wrap(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  // clang-format off
  {"c_cmaes_wrap", (DL_FUNC)&c_cmaes_wrap, 6},
  {NULL, NULL, 0}
  // clang-format on
};

void R_init_libcmaesr(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
