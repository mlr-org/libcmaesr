## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release (new submission NOTE).

## Submission notes

* This is a new release.
* The package bundles a patched copy of the 'libcmaes' C++ library
  (LGPL >= 3). Copyright holders of the bundled code are listed in
  Authors@R with role "cph" and detailed in inst/COPYRIGHTS. The patches
  reroute all library output through Rprintf() and replace C rand() and
  time-based seeding with seeded C++ <random> generators. The vendoring
  script documenting the exact upstream provenance is not part of the
  tarball; it is available at
  https://github.com/mlr-org/libcmaesr/blob/main/tools/vendor-libcmaes.sh
* The method reference for the CMA-ES algorithm (Hansen and Ostermeier,
  2001) is cited in the Description field with a DOI.
* The CRAN status badge in README.md links to
  https://cran.r-project.org/package=libcmaesr, which is currently a 404
  because this is the first submission. The URL resolves once the package
  is accepted.
