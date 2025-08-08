#include <libcmaes/cmaes.h>
#include <Eigen/Dense>
#include "rc_helpers.h"

using namespace libcmaes;
using Eigen::MatrixXd;
using Eigen::VectorXd;


/* 
FIXME:
- allow to set many thing in cmaes params
- map cmaes algo selection to ints, no we can set this as a string
- cmaparams.set_x0(-3.0,3.0); // set x0 domain as [-3,3]^d where d is the dimension of the pro
- look at and test multithreading
- check and implement stop crits
- check what happens if objectives uses wrong data typoes (input and output)
- unit test some NA combos for ctrl params
- show some progess and allow user interupt
- do proper debug printer and put it in rc-helpers
- RC_helpers: add copy function for vecs and matrixes (rows, cols)
- check restarts for ipop and bipop with some form of test 

- provide readme with install instrauctions and minimal example
- do speed test vs 1-2 other packages and maybe perf test
- run R cmd check
- read all files
- read python version

- allow non-vectorized objective. to use multithreading.
- consider noisy case
- look at surrogates
- support unbounded problems


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

  // get dim, lambda, sigma, seed
  // for lambda and sigma, if NA, use default values
  int dim = Rf_length(s_x0);
  int lambda = Rf_asInteger(RC_list_get_el_by_name(s_ctrl, "lambda")); 
  if (lambda == NA_INTEGER) lambda = -1;
  double sigma = Rf_asReal(RC_list_get_el_by_name(s_ctrl, "sigma")); 
  if (R_IsNA(sigma)) sigma = -1;
  int seed = Rf_asInteger(RC_list_get_el_by_name(s_ctrl, "seed")); 
  seed = (seed == NA_INTEGER) ? 0 : seed;
    
 
  // init params and geno-pheno transform, with lin-scaling strategy
  GenoPheno<pwqBoundStrategy, linScalingStrategy> gp(REAL(s_lower), REAL(s_upper), dim);
  CMAParameters<GenoPheno<pwqBoundStrategy,linScalingStrategy>> cmaparams(dim, REAL(s_x0), sigma, lambda, seed, gp);

  // set further params  
  int max_fevals = Rf_asInteger(RC_list_get_el_by_name(s_ctrl, "max_fevals"));
  if (max_fevals != NA_INTEGER) cmaparams.set_max_fevals(max_fevals); 
  int max_iter = Rf_asInteger(RC_list_get_el_by_name(s_ctrl, "max_iter"));
  if (max_iter != NA_INTEGER) cmaparams.set_max_iter(max_iter); 
  double ftarget = Rf_asReal(RC_list_get_el_by_name(s_ctrl, "ftarget"));
  if (!R_IsNA(ftarget)) cmaparams.set_ftarget(ftarget);
  int max_restarts = Rf_asInteger(RC_list_get_el_by_name(s_ctrl, "max_restarts"));
  if (max_restarts != NA_INTEGER) cmaparams.set_restarts(max_restarts);
 
  // dummy needed for constructer, we bypass it later
  FitFunc dummy_scalar = [](const double*, int){ return 0.0; };
  auto es = CMAStrategy<CovarianceUpdate, GenoPheno<pwqBoundStrategy, linScalingStrategy>>
    (dummy_scalar, cmaparams);
  lambda = cmaparams.lambda();  // get lambda from cmaparams if we used -1 for defaults before
  Rprintf("seed: %d; lambda: %d; sigma: %f\n", seed, lambda, sigma);

  SEXP s_x = RC_dblmat_create_PROTECT(lambda, dim); 

  while (!es.stop()) {
    dMat cands_x = es.ask();  // dim x lambda

    // Rprintf("candidates (%ld x %ld):\n", cands_x.rows(), cands_x.cols());
    for (int r=0; r<lambda; r++) {
      dVec x = cands_x.col(r);
      // Rprintf("x(%d): ", r);
      // for (int i = 0; i < dim; i++) {
      //   Rprintf("%f ", x(i));
      // }
      // Rprintf("\n");
      // dVec xp = gp.pheno(x);
      // Rprintf("x-pheno(%d): ", r);
      // for (int i = 0; i < dim; i++) {
      //   Rprintf("%f ", xp(i));
      // }
      // Rprintf("\n");
    }

    // trafo to pheno (!), copy to s_x and eval in R 
    for (int r = 0; r < lambda; r++) {
      dVec x = cands_x.col(r);
      dVec xp = gp.pheno(x);
      for (int i = 0; i < dim; i++) {
        REAL(s_x)[i*lambda + r] = xp(i);
      }
    }
    // Rprintf("s_x copied\n");
    SEXP s_y = RC_tryeval_PROTECT(s_obj, s_x, "libcmaesr: objective evaluation failed!", 1); // on_err: unprotect s_x
    CMASolutions &sols = es.get_solutions();
    for (int r=0; r<sols.size(); r++) {
      sols.get_candidate(r).set_x(cands_x.col(r));
      sols.get_candidate(r).set_fvalue(REAL(s_y)[r]);
    }
    es.update_fevals(sols.size());
    es.tell();     
    es.inc_iter(); 
    UNPROTECT(1); // s_y
  }
 
  // get best seen candidate and trafo to pheno (!)
  const auto &sols = es.get_solutions();
  Candidate bcand = sols.get_best_seen_candidate();
  dVec best_x = gp.pheno(bcand.get_x_dvec()); 

  // copy results to R
  const char* res_names[] = {"x", "y", "edm", "time", "status"};
  SEXP s_res = RC_list_create_withnames_PROTECT(5, res_names);
  SEXP s_res_x = RC_dblvec_create_init_PROTECT(dim, best_x.data());
  SET_VECTOR_ELT(s_res, 0, s_res_x);
  RC_list_set_el_dblscalar(s_res, 1, bcand.get_fvalue());
  RC_list_set_el_dblscalar(s_res, 2, sols.edm());
  RC_list_set_el_dblscalar(s_res, 3, sols.elapsed_time() / 1000.0);
  RC_list_set_el_intscalar(s_res, 4, sols.run_status());
  UNPROTECT(3); // s_x, s_res, s_res_x
    
  return s_res;
}