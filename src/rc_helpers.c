#include "rc_helpers.h"
#include "Rinternals.h"
#include <string.h>

// ********** general **********

void RC_set_names(SEXP s_obj, R_xlen_t n, const char **names) {
  SEXP s_names = PROTECT(Rf_allocVector(STRSXP, n));
  for (R_xlen_t i = 0; i < n; i++) {
    SET_STRING_ELT(s_names, i, Rf_mkChar(names[i]));
  }
  Rf_setAttrib(s_obj, R_NamesSymbol, s_names);
  UNPROTECT(1); // s_names
}

R_xlen_t RC_find_name(SEXP s_obj, const char *name) {
  SEXP s_names = Rf_getAttrib(s_obj, R_NamesSymbol);
  if (!Rf_isString(s_names)) return -1;
  for (R_xlen_t i = 0; i < Rf_xlength(s_names); i++) {
    SEXP s_name = STRING_ELT(s_names, i);
    if (s_name == NA_STRING) continue;
    // normalize to UTF-8 to avoid locale/encoding mismatches
    if (strcmp(Rf_translateCharUTF8(s_name), name) == 0) {
      return i;
    }
  }
  return -1;
}

// ********** scalars **********

const char *RC_charscalar_as_string(SEXP s_x) { return CHAR(STRING_ELT(s_x, 0)); }

// ********** vectors **********

SEXP RC_dblvec_create_init_PROTECT(R_xlen_t n, const double *values) {
  SEXP s_res = PROTECT(Rf_allocVector(REALSXP, n));
  memcpy(REAL(s_res), values, (size_t)n * sizeof(double));
  return s_res;
}

// ********** matrices **********

SEXP RC_dblmat_create_PROTECT(R_xlen_t n_rows, R_xlen_t n_cols) {
  return PROTECT(Rf_allocMatrix(REALSXP, n_rows, n_cols));
}

SEXP RC_dblmat_create_init_PROTECT(R_xlen_t n_rows, R_xlen_t n_cols, const double *x) {
  SEXP s_x = PROTECT(Rf_allocMatrix(REALSXP, n_rows, n_cols));
  memcpy(REAL(s_x), x, n_rows * n_cols * sizeof(double));
  return s_x;
}

// ********** list **********

SEXP RC_list_create_withnames_PROTECT(R_xlen_t n, const char **names) {
  SEXP s_res = PROTECT(Rf_allocVector(VECSXP, n));
  RC_set_names(s_res, n, names);
  return s_res;
}

SEXP RC_list_get_el_by_name(SEXP s_list, const char *name) {
  R_xlen_t i = RC_find_name(s_list, name);
  if (i == -1) return R_NilValue;
  return VECTOR_ELT(s_list, i);
}

SEXP RC_list_set_el_intscalar(SEXP s_list, R_xlen_t idx, r_int32_t x) {
  SET_VECTOR_ELT(s_list, idx, Rf_ScalarInteger(x)); // dont need to protect here
  return s_list;
}

SEXP RC_list_set_el_dblscalar(SEXP s_list, R_xlen_t idx, double x) {
  SET_VECTOR_ELT(s_list, idx, Rf_ScalarReal(x)); // dont need to protect here
  return s_list;
}

SEXP RC_list_set_el_string(SEXP s_list, R_xlen_t idx, const char *x) {
  SET_VECTOR_ELT(s_list, idx, Rf_mkString(x)); // dont need to protect here
  return s_list;
}

// ********** R function calls **********

SEXP RC_tryeval_nothrow_PROTECT(SEXP s_fun, SEXP s_arg, r_int32_t *err) {
  SEXP s_call = PROTECT(Rf_lang2(s_fun, s_arg));
  // returns R_NilValue if there is an error, this is still PROTECTed
  SEXP s_res = PROTECT(R_tryEval(s_call, R_GlobalEnv, err));
  UNPROTECT(1); // s_call
  return s_res;
}

static void RC__chk(void *p) { R_CheckUserInterrupt(); }

r_int32_t RC_interrupt_pending(void) { return !R_ToplevelExec(RC__chk, NULL); }
