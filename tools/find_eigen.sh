#!/usr/bin/env sh
set -e

mkdir -p src

EIGEN_CPPFLAGS=""

if type pkg-config; then
  EIGEN_CPPFLAGS="$(pkg-config --cflags eigen3 2>/dev/null)"
  if [ $? -ne 0 ] && type cmake; then
    EIGEN_CPPFLAGS="$(cmake --find-package -DCOMPILER_ID=GNU -DNAME=Eigen3 -DLANGUAGE=CXX -DMODE=COMPILE 2>/dev/null)"
    if [ $? -ne 0 ]; then
      EIGEN_CPPFLAGS=
    fi
  fi
fi

if [ -z "$EIGEN_CPPFLAGS" ]; then
  cat >&2 <<EOF
ERROR: Eigen3 not found via pkg-config or CMake.
Install Eigen and ensure it is discoverable (PKG_CONFIG_PATH or CMAKE_PREFIX_PATH).
EOF
  exit 1
fi

printf 'EIGEN_CPPFLAGS = %s\n' "$EIGEN_CPPFLAGS" > src/Makevars.eigen
echo "Using system Eigen: $EIGEN_CPPFLAGS"
