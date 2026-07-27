#!/usr/bin/env sh
# Vendor the libcmaes sources into src/libcmaes/.
#
# The sources come from the mlr-org fork of libcmaes, whose main branch
# carries R-specific patches on top of upstream
# https://github.com/CMA-ES/libcmaes:
#   - logging rerouted to Rprintf
#   - IPOP/BIPOP restarts capped by remaining max_fevals
#   - C rand() and unseeded/time-seeded RNGs replaced with seed-derived
#     mt19937 (CRAN compliance, seed reproducibility)
#
# Any further R-specific patches must be pushed to the fork's main branch
# first, then re-vendored here with this script; the vendored copy must stay
# identical to the fork (patches are marked with "libcmaesr patch" comments).
#
# Surrogate model support (surrogatestrategy.*, surrcmaes.h, surrogates/) is
# excluded: it is unused by this package and contains std::cout logging that
# CRAN policy forbids.
#
# The CMake-generated headers cmaes_export.h and libcmaes_config.h are NOT
# copied; static hand-written replacements live in
# src/libcmaes/include/libcmaes/ and must be kept.

set -e

REPO="https://github.com/mlr-org/libcmaes.git"
REF="d86c808cec7088fe6008fc1c57074134ee97b402" # main

DEST="src/libcmaes"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

git clone --quiet "$REPO" "$TMP/libcmaes"
git -C "$TMP/libcmaes" checkout --quiet "$REF"

SOURCES="acovarianceupdate.cc bipopcmastrategy.cc cmaparameters.cc \
  cmasolutions.cc cmastopcriteria.cc cmastrategy.cc covarianceupdate.cc \
  errstats.cc esostrategy.cc ipopcmastrategy.cc pwq_bound_strategy.cc \
  vdcmaupdate.cc"

HEADERS="acovarianceupdate.h bipopcmastrategy.h candidate.h cmaes.h \
  cmaparameters.h cmasolutions.h cmastopcriteria.h cmastrategy.h contour.h \
  covarianceupdate.h eigenmvn.h eo_matrix.h errstats.h esoptimizer.h \
  esostrategy.h genopheno.h ipopcmastrategy.h llogging.h noboundstrategy.h \
  opti_err.h parameters.h pli.h pwq_bound_strategy.h scaling.h vdcmaupdate.h"

mkdir -p "$DEST/include/libcmaes"
for f in $SOURCES; do
  cp "$TMP/libcmaes/src/$f" "$DEST/"
done
for f in $HEADERS; do
  cp "$TMP/libcmaes/include/libcmaes/$f" "$DEST/include/libcmaes/"
done

echo "Vendored libcmaes @ $REF into $DEST"
