#include <atomic>
#include <iostream>
#include <string>
#include <vector>

#include <libcmaes/cmaes.h>
#include <libcmaes/cmasolutions.h>

using namespace libcmaes;

struct EvalCounter {
  std::atomic<long long> count{0};
  double operator()(const double *x, const int &n) {
    // simple sphere; increment true call counter
    count.fetch_add(1, std::memory_order_relaxed);
    double s = 0.0;
    for (int i = 0; i < n; ++i)
      s += x[i] * x[i];
    return s;
  }
};

int main() {
  const int dim = 5;
  std::vector<double> x0(dim, 0.5);
  const double sigma = 0.3;

  // Algorithms to test (restart strategies)
  const std::vector<std::string> algos = {"ipop", "bipop", "sepipop", "sepabipop"};

  for (const std::string &alg : algos) {
    EvalCounter counter;

    CMAParameters<> params(x0, sigma);
    params.set_str_algo(alg);
    params.set_seed(42);
    params.set_quiet(true);
    params.set_mt_feval(false); // ensure single-threaded feval for deterministic counting

    // Give enough budget to allow multiple restarts for restart strategies
    const int max_fevals = 4000;
    params.set_max_fevals(max_fevals);
    params.set_restarts(10);

    FitFunc func = [&](const double *x, const int &n) -> double { return counter(x, n); };

    CMASolutions sols = cmaes<GenoPheno<NoBoundStrategy>>(func, params);

    const long long true_evals = counter.count.load(std::memory_order_relaxed);
    const int sols_fevals = sols.fevals();
    const int sols_nevals = sols.nevals();

    std::cout << "algo=" << alg << ", true_evals=" << true_evals << ", sols.fevals()=" << sols_fevals
              << ", sols.nevals()=" << sols_nevals;

    if (true_evals != sols_fevals || true_evals != sols_nevals) {
      std::cout << " -> MISMATCH";
    }
    std::cout << std::endl;
  }

  return 0;
}
