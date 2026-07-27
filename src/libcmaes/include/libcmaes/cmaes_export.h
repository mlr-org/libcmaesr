/* Static replacement for the CMake-generated export header
 * (GenerateExportHeader). libcmaesr always builds libcmaes as part of the
 * package shared object, so no symbols need to be exported/imported. */

#ifndef CMAES_EXPORT_H
#define CMAES_EXPORT_H

#define CMAES_EXPORT
#define CMAES_NO_EXPORT

#ifndef CMAES_DEPRECATED
#  define CMAES_DEPRECATED
#endif

#ifndef CMAES_DEPRECATED_EXPORT
#  define CMAES_DEPRECATED_EXPORT CMAES_EXPORT CMAES_DEPRECATED
#endif

#ifndef CMAES_DEPRECATED_NO_EXPORT
#  define CMAES_DEPRECATED_NO_EXPORT CMAES_NO_EXPORT CMAES_DEPRECATED
#endif

#endif /* CMAES_EXPORT_H */
