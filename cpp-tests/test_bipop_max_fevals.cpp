#include <iostream>
#include <vector>
#include <cmath>

#include <libcmaes/cmaes.h>
#include <libcmaes/cmasolutions.h>

using namespace libcmaes;

int main() {
  const int dim = 2;
  std::vector<double> x0(dim, 0.5);
  double sigma = 0.2;

  CMAParameters<> params(x0, sigma);
  params.set_str_algo("bipop");
  params.set_seed(42);
  params.set_quiet(false);

  const int max_fevals = 40;
  params.set_max_fevals(max_fevals);
  params.set_restarts(5);


  FitFunc func = [&](const double *x, const int &n) -> double {
    double s = 0.0;
    for (int i = 0; i < n; ++i) s += x[i] * x[i];
    return s;
  };

  CMASolutions sols = cmaes<GenoPheno<NoBoundStrategy>>(func, params);

  return 0;
}
