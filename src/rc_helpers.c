#include "rc_helpers.h"


// ********** general **********

void RC_set_class(SEXP s_obj, const char* class_name) {
  SEXP s_attr = PROTECT(Rf_allocVector(STRSXP, 1));
  SET_STRING_ELT(s_attr, 0, Rf_mkChar(class_name));
  Rf_setAttrib(s_obj, R_ClassSymbol, s_attr);
  UNPROTECT(1); // s_attr
} 

void RC_set_names(SEXP s_obj, int n, const char** names) {
  SEXP s_names = PROTECT(Rf_allocVector(STRSXP, n));
  for (int i = 0; i < n; i++) {
      SET_STRING_ELT(s_names, i, Rf_mkChar(names[i]));
  }
  Rf_setAttrib(s_obj, R_NamesSymbol, s_names);
  UNPROTECT(1); // s_names
}

int RC_find_name(SEXP s_obj, const char *name) {
  SEXP s_names = Rf_getAttrib(s_obj, R_NamesSymbol);
  for (int i = 0; i < Rf_length(s_names); i++) {
    if(strncmp(CHAR(STRING_ELT(s_names, i)), name, strlen(name)) == 0) {
      return i;
    }
  }
  return -1;
}

// ********** scalars **********

const char *RC_charscalar_as_string(SEXP s_x) { 
  return CHAR(STRING_ELT(s_x, 0)); 
}

SEXP RC_intscalar_create_PROTECT(int k) {
  return PROTECT(Rf_ScalarInteger(k));
}

SEXP RC_dblscalar_create_PROTECT(double k) {
  return PROTECT(Rf_ScalarReal(k));
}

// ********** vectors **********

SEXP RC_intvec_create_PROTECT(int n) {
  return PROTECT(Rf_allocVector(INTSXP, n));
}

SEXP RC_dblvec_create_PROTECT(int n) {
  return PROTECT(Rf_allocVector(REALSXP, n));
}

SEXP RC_dblvec_create_init_PROTECT(int n, const double *values) {
  SEXP s_res = PROTECT(Rf_allocVector(REALSXP, n));
  memcpy(REAL(s_res), values, n * sizeof(double));
  return s_res; 
}

// ********** matrices **********

SEXP RC_dblmat_create_PROTECT(int n_rows, int n_cols) {
  return PROTECT(Rf_allocMatrix(REALSXP, n_rows, n_cols));
}

// ********** list **********

SEXP RC_list_create_emptynames_PROTECT(int n) {
  SEXP s_res = PROTECT(Rf_allocVector(VECSXP, n));
  SEXP s_names = PROTECT(Rf_allocVector(STRSXP, n));
  for (int i = 0; i < n; i++) { // initialize names to empty strings
    SET_STRING_ELT(s_names, i, Rf_mkChar(""));
  }
  Rf_setAttrib(s_res, R_NamesSymbol, s_names);
  UNPROTECT(1); // s_names
  return s_res;
}

SEXP RC_list_create_withnames_PROTECT(int n, const char **names) {
  SEXP s_res = PROTECT(Rf_allocVector(VECSXP, n));
  RC_set_names(s_res, n, names);
  return s_res;
}

SEXP RC_list_get_el_by_name(SEXP s_list, const char *name) {
  int i = RC_find_name(s_list, name);
  if (i == -1) return R_NilValue;
  return VECTOR_ELT(s_list, i);
}

// ********** data.frame **********


R_xlen_t RC_df_get_nrows(SEXP s_dt) {
  if (XLENGTH(s_dt) == 0) {
    return 0;
  } else {
    return XLENGTH(VECTOR_ELT(s_dt, 0));
  }
}

SEXP RC_df_get_col_by_name(SEXP s_dt, const char *name) {
  int i = RC_find_name(s_dt, name);
  if (i == -1) return R_NilValue;
  return VECTOR_ELT(s_dt, i);
}

SEXP RC_df_create_allnum_nocolnames_PROTECT(int n_rows, int n_cols) {
  SEXP s_res = PROTECT(Rf_allocVector(VECSXP, n_cols));
  RC_set_class(s_res, "data.frame");
  for (int j = 0; j < n_cols; j++) {
    SET_VECTOR_ELT(s_res, j, Rf_allocVector(REALSXP, n_rows));
  }
  // Set row names to sequential numbers starting from 1
  SEXP s_row_names = PROTECT(Rf_allocVector(INTSXP, n_rows));
  for (int i = 0; i < n_rows; i++) {
      INTEGER(s_row_names)[i] = i + 1;
  }
  Rf_setAttrib(s_res, R_RowNamesSymbol, s_row_names);
  UNPROTECT(1); // s_row_names
  
  return s_res;
}

SEXP RC_df_create_allnum_PROTECT(int n_rows, int n_cols, const char** colnames) {
  SEXP s_res = RC_df_create_allnum_nocolnames_PROTECT(n_rows, n_cols);
  RC_set_names(s_res, n_cols, colnames);
  return s_res;
}


// ********** R6 **********

// this might be conservative and not work for inherited members
SEXP RC_r6_get_member(SEXP s_r6, const char* name) {
  return Rf_findVarInFrame(Rf_install(name), s_r6);
}

void RC_r6_set_member(SEXP s_r6, const char* name, SEXP s_value) {
  Rf_defineVar(Rf_install(name), s_value, s_r6);  // dont need to protect symbols
}

// ********** R function calls **********

SEXP RC_tryeval_PROTECT(SEXP s_fun, SEXP s_arg, const char *errmsg, int on_err_unprotect) {
  SEXP s_call = PROTECT(Rf_lang2(s_fun, s_arg));
  int err = 0;
  SEXP s_res = PROTECT(R_tryEval(s_call, R_GlobalEnv, &err));
  if (err != 0) {
    UNPROTECT(2 + on_err_unprotect); // s_call, s_res, + requested by user
    Rf_error("%s", errmsg);
  }
  UNPROTECT(2); // s_call, s_res
  return s_res;
}
