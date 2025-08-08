#include <libcmaes/cmaes.h>
#include <Eigen/Dense>
#include "rc_helpers.h"

using namespace libcmaes;
using Eigen::MatrixXd;
using Eigen::VectorXd;


/* 
FIXME:
- seed cmaes 
- allow to set many thing in cmaes params
- doc linear scaling
- use lower and upper bounds
  can he have unbounded dec vars? inf bounds?
  for now ONLY support given bounds, but we can also be unbounded later. 
- figure out how we can "disable" option in cmaes params, in R they should then be NA? or NULL?
- map cmaes algo selection to ints, no we can set this as a string
- in docs link to status code of result
- unit test the obj with an exception
- cmaparams.set_x0(-3.0,3.0); // set x0 domain as [-3,3]^d where d is the dimension of the pro
- look at and test multithreading
- cmaes params: set lambda and sigma negative to get default values
- check and implement stop crits
- check via debugging that lin scaling is used
- properly doc the function interface in R. do we even allow a non-vec obj fun?
- create unit tests
- unit test some NA combos for ctrl params
*/

 
//  /* ------------------------------------------------------------------ */
//  /* 2. Factory: map a string to the chosen strategy class               */
//  /* ------------------------------------------------------------------ */
//  template<typename Real = double>
//  std::unique_ptr<ESOStrategy<fit<Real>, model<Real>>>
//  make_strategy(const std::string &variant,
//                FitFunc dummy,               // unused scalar callback
//                CMAParameters<> &params)
//  {
//      using Fit   = fit<Real>;
//      using Model = model<Real>;
//      if (variant == "cma")   return std::make_unique<CMAStrategy      <Fit,Model>>(dummy, params);
//      if (variant == "sep")   return std::make_unique<SepCMAStrategy   <Fit,Model>>(dummy, params);
//      if (variant == "vd")    return std::make_unique<VDCMAStrategy    <Fit,Model>>(dummy, params);
//      if (variant == "ipop")  return std::make_unique<IPOPCMAStrategy  <Fit,Model>>(dummy, params);
//      if (variant == "bipop") return std::make_unique<BIPOPCMAStrategy <Fit,Model>>(dummy, params);
//      throw std::invalid_argument("unknown CMA-ES variant \"" + variant + '"');
//  }
 


extern "C" SEXP c_cmaes_wrap(SEXP s_obj, SEXP s_x0, SEXP s_lower, SEXP s_upper, SEXP s_ctrl) {

  int dim = Rf_length(s_x0);
  double sigma = -1;
  int lambda = -1;
  int seed = Rf_asInteger(RC_list_get_el_by_name(s_ctrl, "seed")); 
  seed = (seed == NA_INTEGER) ? 0 : seed;
  
  
  Rprintf("dim: %d; sigma: %f\n", dim, sigma);
 
  // init params and geno-pheno transform, with lin-scaling strategy
  GenoPheno<pwqBoundStrategy, linScalingStrategy> gp(REAL(s_lower), REAL(s_upper), dim);
  CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>> cmaparams(dim, REAL(s_x0), sigma, lambda, seed, gp);

  // set further params  
  int max_fevals = Rf_asInteger(RC_list_get_el_by_name(s_ctrl, "max_fevals"));
  if (max_fevals != NA_INTEGER) cmaparams.set_max_fevals(max_fevals); 
  int max_iter = Rf_asInteger(RC_list_get_el_by_name(s_ctrl, "max_iter"));
  if (max_iter != NA_INTEGER) cmaparams.set_max_iter(max_iter); 
  double ftarget = Rf_asReal(RC_list_get_el_by_name(s_ctrl, "ftarget"));
  if (ftarget != NA_REAL) cmaparams.set_ftarget(ftarget);
 
  FitFunc dummy_scalar = [](const double*, int){ return 0.0; };
  auto es = CMAStrategy<CovarianceUpdate, GenoPheno<pwqBoundStrategy, linScalingStrategy>>
    (dummy_scalar, cmaparams);
  lambda = cmaparams.lambda();  
  Rprintf("lambda: %d\n", lambda);

  SEXP s_x = RC_dblmat_create_PROTECT(lambda, dim); 

  while (!es.stop()) {
      dMat cands_x = es.ask();  // dim x lambda
      Rprintf("candidates (%ld x %ld):\n", cands_x.rows(), cands_x.cols());
      for(int i = 0; i < cands_x.rows(); i++) {
        for(int j = 0; j < cands_x.cols(); j++) {
          Rprintf("%f ", cands_x(i,j));
        }
        Rprintf("\n");
      }

      // copy data to s_x and eval in R 
      for (int i = 0; i < dim; i++) {
        for (int j = 0; j < lambda; j++) {
          REAL(s_x)[i*lambda + j] = cands_x(i, j);
        }
      }
      Rprintf("s_x copied\n");
      SEXP s_call = PROTECT(Rf_lang2(s_obj, s_x));
      int err = 0;
      SEXP s_y = PROTECT(R_tryEval(s_call, R_GlobalEnv, &err));
      if (err) {
        UNPROTECT(3); // s_x, s_y, s_call
        Rf_error("objective evaluation failed");
      }
      UNPROTECT(2); // s_call, s_y
      // set x and y of candidates
      CMASolutions &sols = es.get_solutions();
      for (int r=0; r<sols.size(); r++) {
  	    sols.get_candidate(r).set_x(cands_x.col(r));
        sols.get_candidate(r).set_fvalue(REAL(s_y)[r]);
      }
      es.update_fevals(sols.size());
      es.tell();     
      es.inc_iter();    
    }
 
  const auto &sols = es.get_solutions();
  Candidate bcand = sols.get_best_seen_candidate();
  dVec best_x = gp.pheno(bcand.get_x_dvec()); 

  // copy results to R
  SEXP s_res_x = PROTECT(Rf_allocVector(REALSXP, dim));
  SEXP s_res_y = PROTECT(Rf_allocVector(REALSXP, 1));
  SEXP s_res_edm = PROTECT(Rf_allocVector(REALSXP, 1));
  SEXP s_res_time = PROTECT(Rf_allocVector(REALSXP, 1));
  SEXP s_res_status = PROTECT(Rf_allocVector(REALSXP, 1));

  for (int i = 0; i < dim; i++) {
      REAL(s_res_x)[i] = best_x[i];
  }
  REAL(s_res_y)[0] = bcand.get_fvalue();
  REAL(s_res_edm)[0] = sols.edm(); 
  REAL(s_res_time)[0] = sols.elapsed_time() / 1000.0; // seconds
  REAL(s_res_status)[0] = sols.run_status(); 

  const char* res_names[] = {"x", "y", "edm", "time", "status"};
  SEXP s_res = RC_list_create_withnames_PROTECT(5, res_names);
  SET_VECTOR_ELT(s_res, 0, s_res_x);
  SET_VECTOR_ELT(s_res, 1, s_res_y);
  SET_VECTOR_ELT(s_res, 2, s_res_edm);
  SET_VECTOR_ELT(s_res, 3, s_res_time);
  SET_VECTOR_ELT(s_res, 4, s_res_status);
  UNPROTECT(7); // s_x, s_res, s_res_x, s_res_y, s_res_edm, s_res_time, s_res_status
    
  return s_res;
}