#ifndef RC_HELPERS_H
#define RC_HELPERS_H

#ifdef __cplusplus
extern "C" {
#endif

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

// ********** general **********

// set class of arbitrary R object
void RC_set_class(SEXP s_obj, const char* class_name);
// set names of arbitrary R object
void RC_set_names(SEXP s_obj, int n, const char** names);
// find name in R_NamesSymbol of arbitrary R object, 
// return index of name, -1 if not found
int RC_find_name(SEXP s_obj, const char *name);

// ********** scalars **********

// convert a charvec(1) to string
const char *RC_charscalar_as_string(SEXP s_x);
// create integer scalar
SEXP RC_intscalar_create_PROTECT(int k);
// create double scalar
SEXP RC_dblscalar_create_PROTECT(double k);

// ********** vectors **********

// create integer vector
SEXP RC_intvec_create_PROTECT(int n);
// create double vector
SEXP RC_dblvec_create_PROTECT(int n);
// create double vector with initial values
SEXP RC_dblvec_create_init_PROTECT(int n, const double *values);

// ********** matrices **********

// create double matrix
SEXP RC_dblmat_create_PROTECT(int n_rows, int n_cols);

// ********** list **********

// create list, named with empty strings
SEXP RC_list_create_emptynames_PROTECT(int n);
// create list, named with provided names
SEXP RC_list_create_withnames_PROTECT(int n, const char **names);
// extract list element by name, returns R_NilValue if not found
SEXP RC_list_get_el_by_name(SEXP s_list, const char *name);

// // set list elements by index
SEXP RC_list_set_el_intscalar(SEXP s_list, int idx, int x);
SEXP RC_list_set_el_dblscalar(SEXP s_list, int idx, double x);


// ********** data.frame **********

// get number of rows in data.frame
R_xlen_t RC_df_get_nrows(SEXP s_dt);
// extract data.frame column by name, returns R_NilValue if not found
SEXP RC_df_get_col_by_name(SEXP s_dt, const char *name);
// create data.frame with all numeric columns, dont set column names
SEXP RC_df_create_allnum_nocolnames_PROTECT(int n_rows, int n_cols);
// create data.frame with all numeric columns, with provided column names
SEXP RC_df_create_allnum_PROTECT(int n_rows, int n_cols, const char** colnames);

// ********** R6 **********

// get R6 member by name
// will return R_UnboundValue when the symbol if symbol is not found in frame
SEXP RC_r6_get_member(SEXP s_r6, const char* name);
// set R6 member by name
// will throw an error if the symbol is not found in frame
void RC_r6_set_member(SEXP s_r6, const char* name, SEXP s_value);

// ********** R function calls **********

// call R function with argument, return result, in simple R_tryEval
// if error, unprotect nr of objects requested by user and throw error with msg
// NB: dont use this if you need more complex cleanup / memory management on error!
// except for the extra unprotect-on-error we assume nothhing more is needed / auto-cleaned-up
SEXP RC_tryeval_PROTECT(SEXP s_fun, SEXP s_arg, const char* errmsg, int on_err_unprotect);


#ifdef __cplusplus
}
#endif

#endif // RC_HELPERS_H