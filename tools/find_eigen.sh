#!/usr/bin/env sh
set -e

mkdir -p src

have() { command -v "$1" >/dev/null 2>&1; }

EIGEN_CFLAGS=""

if have pkg-config; then
  EIGEN_CFLAGS="$(pkg-config --cflags eigen3 2>/dev/null || true)"
fi

if [ -z "$EIGEN_CFLAGS" ] && have cmake; then
  # compiler ID doesn't matter for include discovery; leave generic
  EIGEN_CFLAGS="$(cmake --find-package -DNAME=Eigen3 -DLANGUAGE=CXX -DMODE=COMPILE 2>/dev/null || true)"
fi

# normalize whitespace
EIGEN_CFLAGS=$(printf '%s\n' "$EIGEN_CFLAGS" | tr '\n' ' ' | sed 's/[[:space:]]\+/ /g; s/^ //; s/ $//')

if [ -z "$EIGEN_CFLAGS" ]; then
  cat >&2 <<EOF
ERROR: Eigen3 not found via pkg-config or CMake.
Install Eigen and ensure it is discoverable (PKG_CONFIG_PATH or CMAKE_PREFIX_PATH).
EOF
  exit 1
fi

printf 'PKG_CPPFLAGS += %s\n' "$EIGEN_CFLAGS" > src/Makevars.eigen
echo "Using system Eigen: $EIGEN_CFLAGS"
