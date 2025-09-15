#!/usr/bin/env sh
set -e

EIGEN_CFLAGS=""

# check with "type" if command pkg-config is available
if type pkg-config; then
  EIGEN_CFLAGS="$(pkg-config --cflags eigen3 2>/dev/null)"
  # if exit status is not 0 AND cmake is available, try to find eigen with cmake
  if [ $? -ne 0 ] && type cmake; then
    # COMPILER_ID seems to be necessary, but this should not really influence things
    # in a non-portable way at this stage as we only ask cmake for the include flags
    EIGEN_CFLAGS="$(cmake --find-package -DCOMPILER_ID=GNU -DNAME=Eigen3 -DLANGUAGE=CXX -DMODE=COMPILE 2>/dev/null)"
    # if we fail, set EIGEN_CFLAGS to empty
    if [ $? -ne 0 ]; then
      EIGEN_CFLAGS=
    fi
  fi
fi

if [ -z "$EIGEN_CFLAGS" ]; then
  cat >&2 <<EOF
ERROR: Eigen3 not found via pkg-config or CMake.
Install Eigen and ensure it is discoverable (PKG_CONFIG_PATH or CMAKE_PREFIX_PATH).
EOF
  exit 1
fi

printf 'PKG_CPPFLAGS += %s\n' "$EIGEN_CFLAGS" > src/Makevars.eigen
echo "Using system Eigen: $EIGEN_CFLAGS"
