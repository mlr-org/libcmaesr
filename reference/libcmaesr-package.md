# libcmaesr: R Interface to 'libcmaes'

A lightweight interface to the 'libcmaes' C++ library for the Covariance
Matrix Adaptation Evolution Strategy (CMA-ES). CMA-ES is a
state-of-the-art evolutionary algorithm for the optimization of
difficult non-linear, non-convex black-box functions, as described in
Hansen and Ostermeier (2001)
[doi:10.1162/106365601750190398](https://doi.org/10.1162/106365601750190398)
. Supports the active, separable, and VD (diagonal plus rank-one
covariance) variants of the algorithm as well as the IPOP (increasing
population size) and BIPOP (bi-population) restart strategies. A patched
copy of 'libcmaes' (LGPL \>= 3) is bundled; see the COPYRIGHTS file for
details.

## See also

Useful links:

- <https://github.com/mlr-org/libcmaesr>

- Report bugs at <https://github.com/mlr-org/libcmaesr/issues>

## Author

**Maintainer**: Marc Becker <marcbecker@posteo.de>
([ORCID](https://orcid.org/0000-0002-8115-0400)) \[copyright holder\]

Authors:

- Marc Becker <marcbecker@posteo.de>
  ([ORCID](https://orcid.org/0000-0002-8115-0400)) \[copyright holder\]

- Bernd Bischl <bernd_bischl@gmx.net>
  ([ORCID](https://orcid.org/0000-0001-6002-6980))

- Martin Binder <mlr.developer@mb706.com>

- Lars Kotthoff <lk223@st-andrews.ac.uk>

Other contributors:

- Emmanuel Benazera (author of the bundled 'libcmaes' library)
  \[contributor, copyright holder\]

- Inria (copyright holder of parts of the bundled 'libcmaes' library)
  \[copyright holder\]
