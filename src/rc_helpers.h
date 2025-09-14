#ifndef RC_HELPERS_H
#define RC_HELPERS_H

#ifdef __cplusplus
extern "C" {
#endif

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

// ********** debug prints **********

// Debug printer system - can be switched on/off
#define DEBUG_ENABLED 0 // Set to 1 to enable debug output

#if DEBUG_ENABLED
#define DEBUG_PRINT(fmt, ...) Rprintf(fmt, ##__VA_ARGS__)
#else
#define DEBUG_PRINT(fmt, ...)                                                                                          \
  do {                                                                                                                 \
  } while (0)
#endif

typedef int r_int32_t; // alias for INTSXP payloads

// ********** general **********

// set class of arbitrary R object
// class_name MUST be a null-terminated string
void RC_set_class(SEXP s_obj, const char *class_name);
// set names of arbitrary R object
// names MUST be an array of null-terminated string
void RC_set_names(SEXP s_obj, R_xlen_t n, const char **names);
// find name in R_NamesSymbol of arbitrary R object,
// return index of name, -1 if not found
// name MUST be a null-terminated string
R_xlen_t RC_find_name(SEXP s_obj, const char *name);

// ********** scalars **********

// convert a charvec(1) to string
const char *RC_charscalar_as_string(SEXP s_x);
// create integer scalar
SEXP RC_intscalar_create_PROTECT(r_int32_t k);
// create double scalar
SEXP RC_dblscalar_create_PROTECT(double k);

// ********** vectors **********

// create SEXP vector
SEXP RC_intvec_create_PROTECT(R_xlen_t n);
SEXP RC_dblvec_create_PROTECT(R_xlen_t n);

// create double vector with initial values
SEXP RC_intvec_create_init_PROTECT(R_xlen_t n, const r_int32_t *values);
SEXP RC_dblvec_create_init_PROTECT(R_xlen_t n, const double *values);

// copy to SEXP
void RC_intvec_copy_to_SEXP(const r_int32_t *x, R_xlen_t n, SEXP s_res);
void RC_dblvec_copy_to_SEXP(const double *x, R_xlen_t n, SEXP s_res);

// ********** matrices **********

// create double matrix
SEXP RC_dblmat_create_PROTECT(R_xlen_t n_rows, R_xlen_t n_cols);
// create double matrix with initial values
// NB: the data in x has to be in COLUMN-MAJOR ORDER
SEXP RC_dblmat_create_init_PROTECT(R_xlen_t n_rows, R_xlen_t n_cols, const double *x);

// ********** list **********

// create list, named with empty strings
SEXP RC_list_create_emptynames_PROTECT(R_xlen_t n);
// create list, named with provided names
SEXP RC_list_create_withnames_PROTECT(R_xlen_t n, const char **names);
// extract list element by name, returns R_NilValue if not found
SEXP RC_list_get_el_by_name(SEXP s_list, const char *name);

// // set list elements by index
SEXP RC_list_set_el_intscalar(SEXP s_list, R_xlen_t idx, r_int32_t x);
SEXP RC_list_set_el_dblscalar(SEXP s_list, R_xlen_t idx, double x);
SEXP RC_list_set_el_string(SEXP s_list, R_xlen_t idx, const char *x);

// ********** data.frame **********

// get number of rows in data.frame
R_xlen_t RC_df_get_nrows(SEXP s_dt);
// extract data.frame column by name, returns R_NilValue if not found
SEXP RC_df_get_col_by_name(SEXP s_dt, const char *name);
// create data.frame with all numeric columns, dont set column names
SEXP RC_df_create_allnum_nocolnames_PROTECT(R_xlen_t n_rows, R_xlen_t n_cols);
// create data.frame with all numeric columns, with provided column names
SEXP RC_df_create_allnum_PROTECT(R_xlen_t n_rows, R_xlen_t n_cols, const char **colnames);

// ********** R6 **********

// get R6 member by name
// will return R_UnboundValue when the symbol is not found in frame
SEXP RC_r6_get_member(SEXP s_r6, const char *name);
// set R6 member by name
// will throw an error if the symbol is not found in frame
void RC_r6_set_member(SEXP s_r6, const char *name, SEXP s_value);

// ********** R function calls **********

// call R function with argument, return result, in simple R_tryEval
// if error, unprotect nr of objects requested by user and throw error with msg
// NB: dont use this if you need more complex cleanup / memory management on error!
// except for the extra unprotect-on-error we assume nothing more is needed / auto-cleaned-up
SEXP RC_tryeval_PROTECT(SEXP s_fun, SEXP s_arg, const char *errmsg, r_int32_t on_err_unprotect);
// Try evaluating without raising an error; sets *err to nonzero on error and returns PROTECTed result
// You likely always want to use this instead of RC_tryeval_PROTECT
// NB: result needs to be UNPROTECTed in any case, if there is an error or not
SEXP RC_tryeval_nothrow_PROTECT(SEXP s_fun, SEXP s_arg, r_int32_t *err);
// Check for user interrupt; withpout longjmp for safe interrupt/error handling from C/C++
// Returns 1 if an interrupt is pending
int RC_interrupt_pending(void);

#ifdef __cplusplus
}
#endif

#endif // RC_HELPERS_H
