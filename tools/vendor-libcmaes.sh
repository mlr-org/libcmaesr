#!/usr/bin/env sh
# Vendor the libcmaes sources into src/libcmaes/.
#
# The sources come from our fork of libcmaes, which carries R-specific patches
# on top of upstream https://github.com/CMA-ES/libcmaes:
#   - logging rerouted to Rprintf (branch r_changes)
#   - IPOP/BIPOP restarts capped by remaining max_fevals (branch feat-bipop-budgets)
#
# Additional R-specific patches are applied directly to the vendored copy and
# marked with "libcmaesr patch" comments; re-apply them when re-vendoring, or
# push them to the fork first.
#
# Surrogate model support (surrogatestrategy.*, surrcmaes.h, surrogates/) is
# excluded: it is unused by this package and contains std::cout logging that
# CRAN policy forbids.
#
# The CMake-generated headers cmaes_export.h and libcmaes_config.h are NOT
# copied; static hand-written replacements live in
# src/libcmaes/include/libcmaes/ and must be kept.

set -e

REPO="https://github.com/berndbischl/libcmaes.git"
REF="77c2dc37df650f06dbc6d1e47b9fe905cca30a95" # tip of r_changes

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
echo "Remember to re-apply the 'libcmaesr patch' changes (see git diff)."
