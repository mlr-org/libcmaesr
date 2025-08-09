#include <libcmaes/cmaes.h>
#include <Eigen/Dense>
#include "rc_helpers.h"

using namespace libcmaes;
using Eigen::MatrixXd;
using Eigen::VectorXd;

using MyGenoPheno = GenoPheno<pwqBoundStrategy, linScalingStrategy>;
using MyCMAParameters = CMAParameters<MyGenoPheno>;


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
- allow setting gradients


*/

// temmplated version of the main ask-tell loop, so we can run all templated CMA-ES variants
template <typename Optimizer>
CMASolutions run_optimizer(Optimizer& optimizer, SEXP s_obj, int lambda, int dim, const MyGenoPheno& gp) {
  SEXP s_x = RC_dblmat_create_PROTECT(lambda, dim);
  while (!optimizer.stop()) {
    dMat cands_x = optimizer.ask();
    for (int r = 0; r < lambda; r++) {
      dVec x = cands_x.col(r);
      dVec xp = gp.pheno(x);
      for (int i = 0; i < dim; i++) {
        REAL(s_x)[i*lambda + r] = xp(i);
      }
    }
    SEXP s_y = RC_tryeval_PROTECT(s_obj, s_x, "libcmaesr: objective evaluation failed!", 1); // on_err: unprotect s_x
    CMASolutions &s = optimizer.get_solutions();
    for (int r=0; r<s.size(); r++) {
      s.get_candidate(r).set_x(cands_x.col(r));
      s.get_candidate(r).set_fvalue(REAL(s_y)[r]);
    }
    optimizer.update_fevals(s.size());
    optimizer.tell();
    optimizer.inc_iter();
    UNPROTECT(1); // s_y
  }
  UNPROTECT(1); // s_x
  return optimizer.get_solutions();
}

// I see unfortunatly no way to do this without a switch statement, but at least we can use a template to run all variants
// inspired by libcmaes/cmaes.h
CMASolutions dispatch_optimizer(MyCMAParameters& cmaparams, SEXP s_obj, int lambda, int dim, const MyGenoPheno& gp) {

  // dummy needed for constructer, we bypass it later
  FitFunc dummy_scalar = [](const double*, int){ return 0.0; };
  Rprintf("algo: %d\n", cmaparams.get_algo());
  switch(cmaparams.get_algo()) {
    case CMAES_DEFAULT: {
      ESOptimizer<CMAStrategy<CovarianceUpdate,MyGenoPheno>,MyCMAParameters,CMASolutions> cmaes(dummy_scalar, cmaparams);
      return run_optimizer(cmaes, s_obj, lambda, dim, gp);
    }
    case IPOP_CMAES: {
      ESOptimizer<IPOPCMAStrategy<CovarianceUpdate,MyGenoPheno>,MyCMAParameters,CMASolutions> cmaes(dummy_scalar, cmaparams);
      return run_optimizer(cmaes, s_obj, lambda, dim, gp);
    }
    case BIPOP_CMAES: {
      ESOptimizer<BIPOPCMAStrategy<CovarianceUpdate,MyGenoPheno>,MyCMAParameters,CMASolutions> cmaes(dummy_scalar, cmaparams);
      return run_optimizer(cmaes, s_obj, lambda, dim, gp);
    }
    case aCMAES: {
      ESOptimizer<CMAStrategy<ACovarianceUpdate,MyGenoPheno>,MyCMAParameters,CMASolutions> cmaes(dummy_scalar, cmaparams);
      return run_optimizer(cmaes, s_obj, lambda, dim, gp);
    }
    case aIPOP_CMAES: {
      ESOptimizer<IPOPCMAStrategy<ACovarianceUpdate,MyGenoPheno>,MyCMAParameters,CMASolutions> cmaes(dummy_scalar, cmaparams);
      return run_optimizer(cmaes, s_obj, lambda, dim, gp);
    }
    case aBIPOP_CMAES: {
      ESOptimizer<BIPOPCMAStrategy<ACovarianceUpdate,MyGenoPheno>,MyCMAParameters,CMASolutions> cmaes(dummy_scalar, cmaparams);
      return run_optimizer(cmaes, s_obj, lambda, dim, gp);
    }
    case sepCMAES: {
      cmaparams.set_sep();
      ESOptimizer<CMAStrategy<CovarianceUpdate,MyGenoPheno>,MyCMAParameters,CMASolutions> cmaes(dummy_scalar, cmaparams);
      return run_optimizer(cmaes, s_obj, lambda, dim, gp);
    }
    case sepIPOP_CMAES: {
      cmaparams.set_sep();
      ESOptimizer<IPOPCMAStrategy<CovarianceUpdate,MyGenoPheno>,MyCMAParameters,CMASolutions> cmaes(dummy_scalar, cmaparams);
      return run_optimizer(cmaes, s_obj, lambda, dim, gp);
    }
    case sepBIPOP_CMAES: {
      cmaparams.set_sep();
      ESOptimizer<BIPOPCMAStrategy<CovarianceUpdate,MyGenoPheno>,MyCMAParameters,CMASolutions> cmaes(dummy_scalar, cmaparams);
      return run_optimizer(cmaes, s_obj, lambda, dim, gp);
    }
    case sepaCMAES: {
      cmaparams.set_sep();
      ESOptimizer<CMAStrategy<ACovarianceUpdate,MyGenoPheno>,MyCMAParameters,CMASolutions> cmaes(dummy_scalar, cmaparams);
      return run_optimizer(cmaes, s_obj, lambda, dim, gp);
    }
    case sepaIPOP_CMAES: {
      cmaparams.set_sep();
      ESOptimizer<IPOPCMAStrategy<ACovarianceUpdate,MyGenoPheno>,MyCMAParameters,CMASolutions> cmaes(dummy_scalar, cmaparams);
      return run_optimizer(cmaes, s_obj, lambda, dim, gp);
    }
    case sepaBIPOP_CMAES: {
      cmaparams.set_sep();
      ESOptimizer<BIPOPCMAStrategy<ACovarianceUpdate,MyGenoPheno>,MyCMAParameters,CMASolutions> cmaes(dummy_scalar, cmaparams);
      return run_optimizer(cmaes, s_obj, lambda, dim, gp);
    }
    case VD_CMAES: {
      cmaparams.set_vd();
      ESOptimizer<CMAStrategy<VDCMAUpdate,MyGenoPheno>,MyCMAParameters> cmaes(dummy_scalar, cmaparams);
      return run_optimizer(cmaes, s_obj, lambda, dim, gp);
    }
    case VD_IPOP_CMAES: {
      cmaparams.set_vd();
      ESOptimizer<IPOPCMAStrategy<VDCMAUpdate,MyGenoPheno>,MyCMAParameters> cmaes(dummy_scalar, cmaparams);
      return run_optimizer(cmaes, s_obj, lambda, dim, gp);
    }
    case VD_BIPOP_CMAES: {
      cmaparams.set_vd();
      ESOptimizer<BIPOPCMAStrategy<VDCMAUpdate,MyGenoPheno>,MyCMAParameters> cmaes(dummy_scalar, cmaparams);
      return run_optimizer(cmaes, s_obj, lambda, dim, gp);
    }
    default:
      Rf_error("Unknown CMA-ES variant specified.");
  }
}


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
  MyGenoPheno gp(REAL(s_lower), REAL(s_upper), dim);
  MyCMAParameters cmaparams(dim, REAL(s_x0), sigma, lambda, seed, gp);
  cmaparams.set_str_algo(RC_charscalar_as_string(RC_list_get_el_by_name(s_ctrl, "algo")));

  // set further params  
    int max_fevals = Rf_asInteger(RC_list_get_el_by_name(s_ctrl, "max_fevals"));
  if (max_fevals != NA_INTEGER) cmaparams.set_max_fevals(max_fevals); 
  int max_iter = Rf_asInteger(RC_list_get_el_by_name(s_ctrl, "max_iter"));
  if (max_iter != NA_INTEGER) cmaparams.set_max_iter(max_iter); 
  double ftarget = Rf_asReal(RC_list_get_el_by_name(s_ctrl, "ftarget"));
  if (!R_IsNA(ftarget)) cmaparams.set_ftarget(ftarget);
  int max_restarts = Rf_asInteger(RC_list_get_el_by_name(s_ctrl, "max_restarts"));
  if (max_restarts != NA_INTEGER) cmaparams.set_restarts(max_restarts);
  int elitism = Rf_asInteger(RC_list_get_el_by_name(s_ctrl, "elitism"));
  if (elitism != NA_INTEGER) cmaparams.set_elitism(elitism);
  int tpa = Rf_asInteger(RC_list_get_el_by_name(s_ctrl, "tpa"));
  if (tpa != NA_INTEGER) cmaparams.set_tpa(tpa);
  double dsigma = Rf_asReal(RC_list_get_el_by_name(s_ctrl, "dsigma"));
  if (!R_IsNA(dsigma)) cmaparams.set_tpa_dsigma(dsigma);  

  lambda = cmaparams.lambda();  // get lambda from cmaparams if we used -1 for defaults before
  Rprintf("seed: %d; lambda: %d; sigma: %f\n", seed, lambda, sigma);

  CMASolutions sols = dispatch_optimizer(cmaparams, s_obj, lambda, dim, gp);
 
  // get best seen candidate and trafo to pheno (!)
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
  UNPROTECT(2); // s_res, s_res_x
    
  return s_res;
}