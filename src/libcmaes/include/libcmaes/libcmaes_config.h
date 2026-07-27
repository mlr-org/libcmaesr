/* Static replacement for the CMake-generated configuration header.
 * Only the macros actually referenced by the vendored libcmaes sources are
 * kept: HAVE_DEBUG, HAVE_GLOG and HAVE_SURROG stay undefined (debug logging,
 * Google glog and surrogate models are not used in the R package build). */

#ifndef LIBCMAES_CONFIG_H
#define LIBCMAES_CONFIG_H

#define PACKAGE "libcmaes"
#define PACKAGE_NAME "libcmaes"
#define PACKAGE_STRING "libcmaes 0.10.2"
#define PACKAGE_TARNAME "libcmaes"
#define PACKAGE_VERSION "0.10.2"
#define VERSION "0.10.2"

#endif /* LIBCMAES_CONFIG_H */
