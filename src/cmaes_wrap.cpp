#include "rc_helpers.h"
#include <Eigen/Dense>
#include <libcmaes/cmaes.h>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

using namespace libcmaes;
using Eigen::MatrixXd;
using Eigen::VectorXd;

using MyGenoPheno = GenoPheno<pwqBoundStrategy, linScalingStrategy>;
using MyCMAParameters = CMAParameters<MyGenoPheno>;

// cache for batch evaluation mapping phenotype column pointers to fvalues
static std::unordered_map<const double *, double> G_EVAL_CACHE;
static SEXP G_OBJ;

/*
FIXME:

- in clang-format exlude the libcmaes src dirs

- should be able to reactivate the unizt tests for ipop bibop with budget
problems

- if i runt the unit tests, results dont seem to be completely deterministic
although i seed them? that is a problem..... i could see this here ⠦ | 747 |
single [1] "dim=3, lambda=NA, algo=vdbipopcma" [1] -0.001244280  0.010468126
0.005541255 [1] 0.0001418354 ✖ | 1      871 | single [15.4s]

- use stdint
also this also in RC helpers


- we might have to remnove this from makevars, not portable and warnings in R
CMD check PKG_CXXFLAGS += -Wno-unknown-pragmas

- open issue to overwrite logger better / no exit

- we currently have a fork of libcmaes with 2 branches:
  - r-changes: overwites cout logging to Rprintf
  - feat-bipop-budgets: budget fixes for bipop and ipop

- look at and test multithreading
  openmp i have currently disabled in the build system,
  in any case, we would not be allowed to call into the R Api
  In libcmaes: CMAParameters::set_mt_feval(true/false) toggles OpenMP for
  the population evaluation and finite-difference gradient loops specifically.
  Other OpenMP regions (e.g., in genopheno.h and surrogate rankingsvm.hpp)
  don’t have a library-level switch, so use OMP_NUM_THREADS to control them.

 - we need to ensure that destructors are properly called on exit / exceptoion /
imzterrupt

- at least bipop seems to not respect max_fevals, opened an issue
- provide readme with install instrauctions and minimal example
- do speed test vs 1-2 other packages and maybe perf test
- run R cmd check
- read all files
- read python version

- consider noisy case
  --> there is set_uh in cmaparams
  --> we have to be careful with the global cache here!!!
- look at surrogates
- support unbounded problems
- allow setting gradients


── R CMD check results
──────────────────────────────────────────────────────────────────────────────────────────────────
libcmaesr 0.1 ──── Duration: 3m 25.2s

❯ checking compilation flags in Makevars ... WARNING
  Non-portable flags in variable 'PKG_CXXFLAGS':
    -Wno-unknown-pragmas

❯ checking for GNU extensions in Makefiles ... WARNING
  Found the following file(s) containing GNU extensions:
    src/libcmaes/build/Makefile
    src/libcmaes/build/src/Makefile
  Portable Makefiles do not use GNU extensions such as +=, :=, $(shell),
  $(wildcard), ifeq ... endif, .NOTPARALLEL See section ‘Writing portable
  packages’ in the ‘Writing R Extensions’ manual.

❯ checking compilation flags used ... WARNING
  Compilation used the following non-portable flag(s):
    ‘-Wno-unknown-pragmas’ ‘-mno-omit-leaf-frame-pointer’
  including flag(s) suppressing warnings


❯ checking compiled code ... NOTE
  File ‘libcmaesr/libs/libcmaesr.so’:
    Found ‘_ZSt4cerr’, possibly from ‘std::cerr’ (C++)
      Object: ‘libcmaes/build/src/libcmaes.a’
    Found ‘rand’, possibly from ‘rand’ (C)
      Object: ‘libcmaes/build/src/libcmaes.a’

  Compiled code should not call entry points which might terminate R nor
  write to stdout/stderr instead of to the console, nor use Fortran I/O
  nor system RNGs nor [v]sprintf.

  See ‘Writing portable packages’ in the ‘Writing R Extensions’ manual.




*/

enum class libcmaesr_errcode : int {
  ok = 0,
  bad_return = 1,
  eval_failed = 2,
  user_interrupt = 3
};

// Exception carrying code + context
struct libcmaesr_error : std::runtime_error {
  libcmaesr_errcode code;
  explicit libcmaesr_error(libcmaesr_errcode c, const std::string &msg)
      : std::runtime_error(msg), code(c) {}
};

// Local checker used to throw a typed C++ exception without touching helpers
static inline void check_numvec(SEXP s_y, int expected_length) {
  if (!Rf_isNumeric(s_y)) {
    throw libcmaesr_error(libcmaesr_errcode::bad_return,
                          "libcmaesr: objective must return a numeric vector");
  }
  if (Rf_length(s_y) != expected_length) {
    throw libcmaesr_error(
        libcmaesr_errcode::bad_return,
        std::string(
            "libcmaesr: objective must return a numeric vector of length ") +
            std::to_string(expected_length) + ", got " +
            std::to_string(Rf_length(s_y)));
  }
}

SEXP call_obj_with_error_handling_PROTECT(SEXP s_obj, SEXP s_x, int lambda) {
  int err = 0;
  SEXP s_y = RC_tryeval_nothrow_PROTECT(G_OBJ, s_x, &err);
  if (err != 0) {
    throw libcmaesr_error(libcmaesr_errcode::eval_failed,
                          "libcmaesr: objective evaluation failed!");
  }
  check_numvec(s_y, lambda);
  return s_y;
}

// FitFunc that looks up precomputed fvalues from the cache; falls back to
// single-point R eval
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
    SEXP s_y = call_obj_with_error_handling_PROTECT(G_OBJ, s_x, 1);
    double yval = Rf_asReal(s_y);
    UNPROTECT(2); // s_x, s_y
    return yval;
  }
}

// run strategy with a custom EvalFunc that batch-evaluates the population in R
template <typename Strategy>
static CMASolutions run_with_batch_eval(Strategy &strat, SEXP s_obj) {
  EvalFunc evalf = [&](const dMat &cands, const dMat &phenocands) {
    const int dim = phenocands.rows();
    const int lambda = phenocands.cols();

    // allow user interrupt before allocating/protecting
    if (RC_interrupt_pending())
      throw libcmaesr_error(libcmaesr_errcode::user_interrupt,
                            "libcmaesr: user interrupt");

    // create lambda x dim matrix for R objective, fill with phenotype
    SEXP s_x = RC_dblmat_create_PROTECT(lambda, dim);

    // Column-major map over R memory, unaligned for safety
    Eigen::Map<
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>,
        Eigen::Unaligned>
        X(REAL(s_x), lambda, dim);
    // Transpose phenocands (dim x lambda) into X (lambda x dim)
    X = phenocands.transpose();

    // call R once on the whole batch
    SEXP s_y = call_obj_with_error_handling_PROTECT(s_obj, s_x, lambda);

    // populate cache: map phenotype column pointer -> fvalue
    G_EVAL_CACHE.clear();
    // FIXME: should we do this here? just do this once globally?
    G_EVAL_CACHE.reserve(lambda);
    for (int r = 0; r < lambda; ++r) {
      G_EVAL_CACHE.emplace(phenocands.col(r).data(), REAL(s_y)[r]);
    }
    // delegate to internal eval to mirror all behaviors (UH, elitism, fevals,
    // etc.)
    strat.eval(cands, phenocands);
    UNPROTECT(2); // s_x, s_y
  };

  AskFunc askf = [&]() { return strat.ask(); };
  TellFunc tellf = [&]() { strat.tell(); };

  strat.optimize(evalf, askf, tellf);
  return strat.get_solutions();
}

// minimal dispatch that mirrors libcmaes::cmaes but routes through
// run_with_batch_eval
static CMASolutions dispatch_with_batch_eval(MyCMAParameters &cmaparams,
                                             SEXP s_obj) {
  // FitFunc that looks up precomputed fvalues from the cache; falls back to
  // single-point R eval
  FitFunc cachedFF = [](const double *x, const int &n) {
    return cached_fitfunc_impl(x, n);
  };

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

std::pair<MyCMAParameters, MyGenoPheno> cmaes_setup(SEXP s_x0, SEXP s_lower,
                                                    SEXP s_upper, SEXP s_ctrl) {
  int dim = Rf_length(s_x0);
  int lambda = Rf_asInteger(RC_list_get_el_by_name(s_ctrl, "lambda"));
  if (lambda == NA_INTEGER)
    lambda = -1;
  double sigma = Rf_asReal(RC_list_get_el_by_name(s_ctrl, "sigma"));
  if (R_IsNA(sigma))
    sigma = -1;
  int seed = Rf_asInteger(RC_list_get_el_by_name(s_ctrl, "seed"));
  seed = (seed == NA_INTEGER) ? 0 : seed;

  DEBUG_PRINT("cmaes_setup: dim: %d; lambda: %d; sigma: %f; seed: %d\n", dim,
              lambda, sigma, seed);

  // init params and geno-pheno transform, with lin-scaling strategy
  MyGenoPheno gp(REAL(s_lower), REAL(s_upper), dim);
  MyCMAParameters cmaparams(dim, REAL(s_x0), sigma, lambda, seed, gp);

  // set further params
  cmaparams.set_maximize(
      Rf_asInteger(RC_list_get_el_by_name(s_ctrl, "maximize")));
  cmaparams.set_str_algo(
      RC_charscalar_as_string(RC_list_get_el_by_name(s_ctrl, "algo")));
  int max_fevals = Rf_asInteger(RC_list_get_el_by_name(s_ctrl, "max_fevals"));
  if (max_fevals != NA_INTEGER)
    cmaparams.set_max_fevals(max_fevals);
  int max_iter = Rf_asInteger(RC_list_get_el_by_name(s_ctrl, "max_iter"));
  if (max_iter != NA_INTEGER)
    cmaparams.set_max_iter(max_iter);
  double ftarget = Rf_asReal(RC_list_get_el_by_name(s_ctrl, "ftarget"));
  if (!R_IsNA(ftarget))
    cmaparams.set_ftarget(ftarget);
  double f_tolerance = Rf_asReal(RC_list_get_el_by_name(s_ctrl, "f_tolerance"));
  if (!R_IsNA(f_tolerance))
    cmaparams.set_ftolerance(f_tolerance);
  double x_tolerance = Rf_asReal(RC_list_get_el_by_name(s_ctrl, "x_tolerance"));
  if (!R_IsNA(x_tolerance))
    cmaparams.set_xtolerance(x_tolerance);
  int max_restarts =
      Rf_asInteger(RC_list_get_el_by_name(s_ctrl, "max_restarts"));
  if (max_restarts != NA_INTEGER)
    cmaparams.set_restarts(max_restarts);
  int elitism = Rf_asInteger(RC_list_get_el_by_name(s_ctrl, "elitism"));
  if (elitism != NA_INTEGER)
    cmaparams.set_elitism(elitism);
  int tpa = Rf_asInteger(RC_list_get_el_by_name(s_ctrl, "tpa"));
  if (tpa != NA_INTEGER)
    cmaparams.set_tpa(tpa);
  double dsigma = Rf_asReal(RC_list_get_el_by_name(s_ctrl, "tpa_dsigma"));
  if (!R_IsNA(dsigma))
    cmaparams.set_tpa_dsigma(dsigma);
  bool quiet = Rf_asLogical(RC_list_get_el_by_name(s_ctrl, "quiet"));
  cmaparams.set_quiet(quiet);
  SEXP s_x0_lower = RC_list_get_el_by_name(s_ctrl, "x0_lower");
  SEXP s_x0_upper = RC_list_get_el_by_name(s_ctrl, "x0_upper");
  if (s_x0_lower != R_NilValue && s_x0_upper != R_NilValue) {
    cmaparams.set_x0(REAL(s_x0_lower), REAL(s_x0_upper));
  }
  cmaparams.set_mt_feval(
      false); // disable openmp, we are not allowed to call into R api
  DEBUG_PRINT("cmaes_setup: done\n");

  return std::make_pair(cmaparams, gp);
}

extern "C" SEXP c_cmaes_wrap(SEXP s_obj, SEXP s_x0, SEXP s_lower, SEXP s_upper,
                             SEXP s_ctrl, SEXP s_batch) {
  try {
    std::pair<MyCMAParameters, MyGenoPheno> setup =
        cmaes_setup(s_x0, s_lower, s_upper, s_ctrl);
    MyCMAParameters &cmaparams = setup.first;
    MyGenoPheno &gp = setup.second;
    CMASolutions sols;
    bool batch = Rf_asLogical(s_batch);
    G_OBJ = s_obj;
    if (!batch) {
      FitFunc func = [](const double *x, const int &n) -> double {
        // allow user interrupt before allocating/protecting
        if (RC_interrupt_pending())
          throw libcmaesr_error(libcmaesr_errcode::user_interrupt,
                                "libcmaesr: user interrupt");
        // setup R dbl vec, copy, eval, return
        SEXP s_x = RC_dblvec_create_init_PROTECT(n, x);
        SEXP s_y = call_obj_with_error_handling_PROTECT(G_OBJ, s_x, 1);
        double y = Rf_asReal(s_y);
        UNPROTECT(2); // s_x, s_y
        return y;
      };
      // now run via libcmaes::cmaes servive fun, easy
      sols = cmaes<MyGenoPheno>(func, cmaparams);
    } else {
      sols = dispatch_with_batch_eval(cmaparams, s_obj);
    }

    // get best seen candidate and transform to pheno (!)
    Candidate bcand = sols.get_best_seen_candidate();
    dVec best_x = gp.pheno(bcand.get_x_dvec());
    double best_y = bcand.get_fvalue();
    // maximize -> negate best_y which was internally multiplied by -1
    best_y = cmaparams.get_maximize() ? -best_y : best_y;

    // copy results to R
    const char *res_names[] = {"x",    "y",           "edm",
                               "time", "status_code", "status_msg"};
    SEXP s_res = RC_list_create_withnames_PROTECT(6, res_names);
    SEXP s_res_x = RC_dblvec_create_init_PROTECT(best_x.size(), best_x.data());
    SET_VECTOR_ELT(s_res, 0, s_res_x);
    RC_list_set_el_dblscalar(s_res, 1, best_y);
    RC_list_set_el_dblscalar(s_res, 2, sols.edm());
    RC_list_set_el_dblscalar(s_res, 3, sols.elapsed_time() / 1000.0); // in secs
    RC_list_set_el_intscalar(s_res, 4, sols.run_status());
    SET_VECTOR_ELT(s_res, 5, Rf_mkString(sols.status_msg().c_str()));
    UNPROTECT(2); // s_res, s_res_x
    return s_res;
  } catch (const libcmaesr_error &e) {
    Rf_error("%s", e.what());
  } catch (...) {
    Rf_error("libcmaesr: unknown error");
  }
  // not reached, keeps compiler happy
  return R_NilValue;
}
