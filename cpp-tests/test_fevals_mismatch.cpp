#include <iostream>
#include <string>
#include <vector>

#include <libcmaes/cmaes.h>
#include <libcmaes/cmasolutions.h>

using namespace libcmaes;

static int COUNTER = 0;

int main() {
  const int dim = 5;
  std::vector<double> x0(dim, 0.5);
  const double sigma = 0.3;

  // Algorithms to test (restart strategies)
  const std::vector<std::string> algos = {"ipop", "bipop", "sepipop", "sepabipop"};

  for (const std::string &alg : algos) {
    COUNTER = 0;
    CMAParameters<> params(x0, sigma);
    params.set_str_algo(alg);
    params.set_seed(42);
    params.set_mt_feval(false); // ensure single-threaded feval

    // Give enough budget to allow multiple restarts for restart strategies
    params.set_max_fevals(4000);
    params.set_restarts(10);

    FitFunc func = [&](const double *x, const int &n) -> double {
      COUNTER++;
      return 1.0;
    };

    CMASolutions sols = cmaes<GenoPheno<NoBoundStrategy>>(func, params);

    std::cout << "algo=" << alg << ", true_evals=" << COUNTER << ", sols.fevals() = " << sols.fevals()
              << ", sols.nevals() = " << sols.nevals() << std::endl;
  }

  return 0;
}
