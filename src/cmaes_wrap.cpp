#include <asm-generic/errno.h>
#include <libcmaes/cmaes.h>
#include <Eigen/Dense>
#include <unordered_map>
#include "rc_helpers.h"

using namespace libcmaes;
using Eigen::MatrixXd;
using Eigen::VectorXd;

using MyGenoPheno = GenoPheno<pwqBoundStrategy, linScalingStrategy>;
using MyCMAParameters = CMAParameters<MyGenoPheno>;

// cache for batch evaluation mapping phenotype column pointers to fvalues
static std::unordered_map<const double*, double> G_EVAL_CACHE;
static SEXP G_OBJ;

/*
FIXME:
- cmaparams.set_x0(-3.0,3.0); // set x0 domain as [-3,3]^d where d is the dimension of the pro
- look at and test multithreading
- check what happens if objectives uses wrong data typoes (input and output)
- unit test some NA combos for ctrl params
- RC_helpers: add copy function for vecs and matrixes (rows, cols)
- add lic etc for libcmaes
- check if we can make teh cache eval faster
- check restarts for ipop and bipop with some form of test
- reenable restart tests in test_batch.R

- at least bipop seems to not respect max_fevals, opened an issue
- provide readme with install instrauctions and minimal example
- do speed test vs 1-2 other packages and maybe perf test
- run R cmd check
- read all files
- read python version

- allow non-vectorized objective. to use multithreading.
- consider noisy case
  --> there is set_uh in cmaparams
- look at surrogates
- support unbounded problems
- allow setting gradients


*/

// FitFunc that looks up precomputed fvalues from the cache; falls back to single-point R eval
static double cached_fitfunc_impl(const double *x, int dim) {
  DEBUG_PRINT("cached_fitfunc_impl: dim=%d\n", dim);
  // FIXME can we make this faster?
  auto it = G_EVAL_CACHE.find(x);
  if (it != G_EVAL_CACHE.end()) {
    return it->second;
  } else {
    DEBUG_PRINT("cached_fitfunc_impl: not found\n");
    // create matrix with 1 row and eval in R
    SEXP s_x = RC_dblmat_create_init_PROTECT(1, dim, x);
    SEXP s_y = RC_tryeval_PROTECT(G_OBJ, s_x, "libcmaesr: objective evaluation failed!", 1); // unprotect s_x on err
    UNPROTECT(2); // s_x, s_y
    return REAL(s_y)[0];
  }
}

// run strategy with a custom EvalFunc that batch-evaluates the population in R
template <typename Strategy>
static CMASolutions run_with_batch_eval(Strategy &strat, SEXP s_obj) {
  EvalFunc evalf = [&](const dMat &cands, const dMat &phenocands) {
    const int dim = phenocands.rows();
    const int lambda = phenocands.cols();

    // create lambda x dim matrix for R objective, fill with phenotype
    SEXP s_x = RC_dblmat_create_PROTECT(lambda, dim);

    // Column-major map over R memory, unaligned for safety
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>, Eigen::Unaligned>
        X(REAL(s_x), lambda, dim);
    // Transpose phenocands (dim x lambda) into X (lambda x dim)
    X = phenocands.transpose();

    // call R once on the whole batch
    SEXP s_y = RC_tryeval_PROTECT(s_obj, s_x, "libcmaesr: objective evaluation failed!", 1); // unprotect s_x on err

    // populate cache: map phenotype column pointer -> fvalue
    G_EVAL_CACHE.clear();
    for (int r = 0; r < lambda; ++r) {
        G_EVAL_CACHE.emplace(phenocands.col(r).data(), REAL(s_y)[r]);
    }
    // delegate to internal eval to mirror all behaviors (UH, elitism, fevals, etc.)
    strat.eval(cands, phenocands);
    UNPROTECT(2); // s_x, s_y
  };

  AskFunc askf = [&]() { return strat.ask(); };
  TellFunc tellf = [&]() { strat.tell(); };

  strat.optimize(evalf, askf, tellf);
  return strat.get_solutions();
}

// minimal dispatch that mirrors libcmaes::cmaes but routes through run_with_batch_eval
static CMASolutions dispatch_with_batch_eval(MyCMAParameters &cmaparams, SEXP s_obj) {
    // FitFunc that looks up precomputed fvalues from the cache; falls back to single-point R eval
    FitFunc cachedFF = [](const double *x, const int &n) { return cached_fitfunc_impl(x, n); };

  switch (cmaparams.get_algo()) {
    case CMAES_DEFAULT: {
      CMAStrategy<CovarianceUpdate, MyGenoPheno> strat(cachedFF, cmaparams);
      return run_with_batch_eval(strat, s_obj);
    }
    case IPOP_CMAES: {
      IPOPCMAStrategy<CovarianceUpdate, MyGenoPheno> strat(cachedFF, cmaparams);
      return run_with_batch_eval(strat, s_obj);
    }
    case BIPOP_CMAES: {
      BIPOPCMAStrategy<CovarianceUpdate, MyGenoPheno> strat(cachedFF, cmaparams);
      return run_with_batch_eval(strat, s_obj);
    }
    case aCMAES: {
      CMAStrategy<ACovarianceUpdate, MyGenoPheno> strat(cachedFF, cmaparams);
      return run_with_batch_eval(strat, s_obj);
    }
    case aIPOP_CMAES: {
      IPOPCMAStrategy<ACovarianceUpdate, MyGenoPheno> strat(cachedFF, cmaparams);
      return run_with_batch_eval(strat, s_obj);
    }
    case aBIPOP_CMAES: {
      BIPOPCMAStrategy<ACovarianceUpdate, MyGenoPheno> strat(cachedFF, cmaparams);
      return run_with_batch_eval(strat, s_obj);
    }
    case sepCMAES: {
      cmaparams.set_sep();
      CMAStrategy<CovarianceUpdate, MyGenoPheno> strat(cachedFF, cmaparams);
      return run_with_batch_eval(strat, s_obj);
    }
    case sepIPOP_CMAES: {
      cmaparams.set_sep();
      IPOPCMAStrategy<CovarianceUpdate, MyGenoPheno> strat(cachedFF, cmaparams);
      return run_with_batch_eval(strat, s_obj);
    }
    case sepBIPOP_CMAES: {
      cmaparams.set_sep();
      BIPOPCMAStrategy<CovarianceUpdate, MyGenoPheno> strat(cachedFF, cmaparams);
      return run_with_batch_eval(strat, s_obj);
    }
    case sepaCMAES: {
      cmaparams.set_sep();
      CMAStrategy<ACovarianceUpdate, MyGenoPheno> strat(cachedFF, cmaparams);
      return run_with_batch_eval(strat, s_obj);
    }
    case sepaIPOP_CMAES: {
      cmaparams.set_sep();
      IPOPCMAStrategy<ACovarianceUpdate, MyGenoPheno> strat(cachedFF, cmaparams);
      return run_with_batch_eval(strat, s_obj);
    }
    case sepaBIPOP_CMAES: {
      cmaparams.set_sep();
      BIPOPCMAStrategy<ACovarianceUpdate, MyGenoPheno> strat(cachedFF, cmaparams);
      return run_with_batch_eval(strat, s_obj);
    }
    case VD_CMAES: {
      cmaparams.set_vd();
      CMAStrategy<VDCMAUpdate, MyGenoPheno> strat(cachedFF, cmaparams);
      return run_with_batch_eval(strat, s_obj);
    }
    case VD_IPOP_CMAES: {
      cmaparams.set_vd();
      IPOPCMAStrategy<VDCMAUpdate, MyGenoPheno> strat(cachedFF, cmaparams);
      return run_with_batch_eval(strat, s_obj);
    }
    case VD_BIPOP_CMAES: {
      cmaparams.set_vd();
      BIPOPCMAStrategy<VDCMAUpdate, MyGenoPheno> strat(cachedFF, cmaparams);
      return run_with_batch_eval(strat, s_obj);
    }
    default:
      Rf_error("Unknown CMA-ES variant specified.");
  }

  // not reached
  return CMASolutions();
}


std::pair<MyCMAParameters, MyGenoPheno> cmaes_setup(SEXP s_x0, SEXP s_lower, SEXP s_upper, SEXP s_ctrl) {
    int dim = Rf_length(s_x0);
    int lambda = Rf_asInteger(RC_list_get_el_by_name(s_ctrl, "lambda"));
    if (lambda == NA_INTEGER) lambda = -1;
    double sigma = Rf_asReal(RC_list_get_el_by_name(s_ctrl, "sigma"));
    if (R_IsNA(sigma)) sigma = -1;
    int seed = Rf_asInteger(RC_list_get_el_by_name(s_ctrl, "seed"));
    seed = (seed == NA_INTEGER) ? 0 : seed;

    DEBUG_PRINT("cmaes_setup: dim: %d; lambda: %d; sigma: %f; seed: %d\n", dim, lambda, sigma, seed);

    // init params and geno-pheno transform, with lin-scaling strategy
    MyGenoPheno gp(REAL(s_lower), REAL(s_upper), dim);
    MyCMAParameters cmaparams(dim, REAL(s_x0), sigma, lambda, seed, gp);

    // set further params
    cmaparams.set_maximize(Rf_asInteger(RC_list_get_el_by_name(s_ctrl, "maximize")));
    cmaparams.set_str_algo(RC_charscalar_as_string(RC_list_get_el_by_name(s_ctrl, "algo")));
    int max_fevals = Rf_asInteger(RC_list_get_el_by_name(s_ctrl, "max_fevals"));
    if (max_fevals != NA_INTEGER) cmaparams.set_max_fevals(max_fevals);
    int max_iter = Rf_asInteger(RC_list_get_el_by_name(s_ctrl, "max_iter"));
    if (max_iter != NA_INTEGER) cmaparams.set_max_iter(max_iter);
    double ftarget = Rf_asReal(RC_list_get_el_by_name(s_ctrl, "ftarget"));
    if (!R_IsNA(ftarget)) cmaparams.set_ftarget(ftarget);
    double f_tolerance = Rf_asReal(RC_list_get_el_by_name(s_ctrl, "f_tolerance"));
    if (!R_IsNA(f_tolerance)) cmaparams.set_ftolerance(f_tolerance);
    double x_tolerance = Rf_asReal(RC_list_get_el_by_name(s_ctrl, "x_tolerance"));
    if (!R_IsNA(x_tolerance)) cmaparams.set_xtolerance(x_tolerance);
    int max_restarts = Rf_asInteger(RC_list_get_el_by_name(s_ctrl, "max_restarts"));
    if (max_restarts != NA_INTEGER) cmaparams.set_restarts(max_restarts);
    int elitism = Rf_asInteger(RC_list_get_el_by_name(s_ctrl, "elitism"));
    if (elitism != NA_INTEGER) cmaparams.set_elitism(elitism);
    int tpa = Rf_asInteger(RC_list_get_el_by_name(s_ctrl, "tpa"));
    if (tpa != NA_INTEGER) cmaparams.set_tpa(tpa);
    double dsigma = Rf_asReal(RC_list_get_el_by_name(s_ctrl, "tpa_dsigma"));
    if (!R_IsNA(dsigma)) cmaparams.set_tpa_dsigma(dsigma);
    bool quiet = Rf_asLogical(RC_list_get_el_by_name(s_ctrl, "quiet"));
    cmaparams.set_quiet(quiet);

    DEBUG_PRINT("cmaes_setup: done\n");

    return std::make_pair(cmaparams, gp);
}

SEXP create_SEXP_result(CMASolutions& sols, MyGenoPheno& gp, MyCMAParameters& cmaparams) {
    // get best seen candidate and transform to pheno (!)
    Candidate bcand = sols.get_best_seen_candidate();
    dVec best_x = gp.pheno(bcand.get_x_dvec());
    double best_y = bcand.get_fvalue();
    // maximize -> negate best_y which was internally multiplied by -1
    best_y = cmaparams.get_maximize() ? -best_y : best_y;

    // copy results to R
    const char* res_names[] = {"x", "y", "edm", "time", "status"};
    SEXP s_res = RC_list_create_withnames_PROTECT(5, res_names);
    SEXP s_res_x = RC_dblvec_create_init_PROTECT(best_x.size(), best_x.data());
    SET_VECTOR_ELT(s_res, 0, s_res_x);
    RC_list_set_el_dblscalar(s_res, 1, best_y);
    RC_list_set_el_dblscalar(s_res, 2, sols.edm());
    RC_list_set_el_dblscalar(s_res, 3, sols.elapsed_time() / 1000.0); // in secs
    RC_list_set_el_intscalar(s_res, 4, sols.run_status());
    UNPROTECT(2); // s_res, s_res_x
    return s_res;
}


extern "C" SEXP c_cmaes_wrap_single(SEXP s_obj, SEXP s_x0, SEXP s_lower, SEXP s_upper, SEXP s_ctrl) {

  auto [cmaparams, gp] = cmaes_setup(s_x0, s_lower, s_upper, s_ctrl);

  // scalar fitfunc, simple case
  FitFunc func = [&](const double *x, const int &n) -> double {
    // setup R dbl vec, copy, eval, return
    SEXP s_x = RC_dblvec_create_init_PROTECT(n, x);
    SEXP s_y = RC_tryeval_PROTECT(s_obj, s_x, "libcmaesr: objective evaluation failed!", 1); // unprotect s_x on err
    UNPROTECT(2); // s_x, s_y
    return Rf_asReal(s_y);
  };

  // now run via libcmaes::cmaes servive fun, easy
  CMASolutions sols = cmaes<MyGenoPheno>(func, cmaparams);
  // max_mult and manual pheno not needed for single
  return create_SEXP_result(sols, gp, cmaparams);
}



extern "C" SEXP c_cmaes_wrap_batch(SEXP s_obj, SEXP s_x0, SEXP s_lower, SEXP s_upper, SEXP s_ctrl) {
    auto [cmaparams, gp] = cmaes_setup(s_x0, s_lower, s_upper, s_ctrl);
    G_OBJ = s_obj;
    CMASolutions sols = dispatch_with_batch_eval(cmaparams, s_obj);
    return create_SEXP_result(sols, gp, cmaparams);
}
